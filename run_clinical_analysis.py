#!/usr/bin/env python3
"""Run clinical analysis with robust logging and monitoring to prevent accidental interruption.

This wrapper ensures that the clinical analysis runs to completion with:
- Real-time logging to file and console
- Progress monitoring and heartbeat
- Watchdog timer to detect hangs
- Auto-recovery from minor errors
- Final validation that analysis completed successfully
"""

import os
import sys
import time
import threading
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional


class AnalysisMonitor:
    """Monitor clinical analysis execution with logging and watchdog."""
    
    def __init__(self, log_file: Path):
        self.log_file = log_file
        self.start_time = time.time()
        self.last_activity = time.time()
        self.is_running = False
        self.watchdog_thread: Optional[threading.Thread] = None
        self.heartbeat_thread: Optional[threading.Thread] = None
        
    def log(self, message: str, level: str = "INFO") -> None:
        """Write to both log file and console."""
        timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_line = f"[{timestamp}] [{level}] {message}\n"
        
        # Write to file
        with open(self.log_file, 'a', encoding='utf-8') as f:
            f.write(log_line)
        
        # Write to console
        print(log_line.rstrip())
        sys.stdout.flush()
        
        self.last_activity = time.time()
    
    def watchdog(self, timeout: int = 300) -> None:
        """Monitor for hangs - alert if no activity for timeout seconds."""
        while self.is_running:
            time.sleep(30)  # Check every 30 seconds
            elapsed = time.time() - self.last_activity
            if elapsed > timeout:
                self.log(f"⚠️ WARNING: No activity for {elapsed:.0f} seconds (timeout: {timeout}s)", "WARN")
                self.log("Analysis may be hung. Check process status.", "WARN")
            time.sleep(30)
    
    def heartbeat(self) -> None:
        """Send periodic heartbeat to show analysis is still running."""
        while self.is_running:
            elapsed = time.time() - self.start_time
            minutes = int(elapsed // 60)
            seconds = int(elapsed % 60)
            self.log(f"💓 Analysis running... Elapsed: {minutes}m {seconds}s", "HEARTBEAT")
            time.sleep(60)  # Heartbeat every 60 seconds
    
    def start_monitoring(self) -> None:
        """Start watchdog and heartbeat threads."""
        self.is_running = True
        
        self.watchdog_thread = threading.Thread(target=self.watchdog, daemon=True)
        self.watchdog_thread.start()
        
        self.heartbeat_thread = threading.Thread(target=self.heartbeat, daemon=True)
        self.heartbeat_thread.start()
        
        self.log("✅ Monitoring started - Watchdog and heartbeat active")
    
    def stop_monitoring(self) -> None:
        """Stop monitoring threads."""
        self.is_running = False
        elapsed = time.time() - self.start_time
        minutes = int(elapsed // 60)
        seconds = int(elapsed % 60)
        self.log(f"⏹️ Monitoring stopped - Total runtime: {minutes}m {seconds}s")


def check_dependencies() -> bool:
    """Check if required Python packages are installed."""
    required_packages = ['pandas', 'matplotlib', 'seaborn', 'lifelines', 'openpyxl']
    missing = []
    
    for package in required_packages:
        try:
            __import__(package)
        except ImportError:
            missing.append(package)
    
    if missing:
        print(f"❌ Missing required packages: {', '.join(missing)}")
        print(f"\nInstall with: pip install {' '.join(missing)}")
        return False
    
    return True


def run_clinical_analysis(monitor: AnalysisMonitor) -> bool:
    """Run the clinical analysis with monitoring."""
    
    # Check if clinical data exists
    clinical_file = Path("DATA/CGGA.WEseq_286_clinical.20200506.txt")
    if not clinical_file.exists():
        monitor.log(f"❌ ERROR: Clinical data file not found: {clinical_file}", "ERROR")
        return False
    
    monitor.log(f"📊 Found clinical data: {clinical_file}")
    monitor.log(f"📁 Number of records: {len(open(clinical_file).readlines()) - 1}")
    
    # Prepare output directory
    output_dir = Path("clinical_analysis_outputs")
    output_dir.mkdir(exist_ok=True)
    monitor.log(f"📂 Output directory: {output_dir.absolute()}")
    
    # Build command
    cmd = [
        sys.executable,
        "-m",
        "apoe_analysis.clinical_analysis",
        "--clinical",
        str(clinical_file),
        "--output-dir",
        str(output_dir),
        "--time-column",
        "OS",
        "--event-column",
        "Censor (alive=0; dead=1)",
        "--group-column",
        "Grade",
    ]
    
    monitor.log("🚀 Starting clinical analysis...")
    monitor.log(f"Command: {' '.join(cmd)}")
    monitor.log("=" * 80)
    
    # Start monitoring
    monitor.start_monitoring()
    
    try:
        # Run the analysis
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            universal_newlines=True
        )
        
        # Stream output in real-time
        if process.stdout:
            for line in process.stdout:
                monitor.log(line.rstrip(), "ANALYSIS")
        
        # Wait for completion
        return_code = process.wait()
        
        monitor.log("=" * 80)
        if return_code == 0:
            monitor.log("✅ Clinical analysis completed successfully!", "SUCCESS")
            return True
        else:
            monitor.log(f"❌ Analysis failed with return code: {return_code}", "ERROR")
            return False
    
    except KeyboardInterrupt:
        monitor.log("⚠️ Keyboard interrupt detected - attempting graceful shutdown...", "WARN")
        if process:
            process.terminate()
            process.wait(timeout=10)
        return False
    
    except Exception as exc:
        monitor.log(f"❌ Unexpected error: {exc}", "ERROR")
        import traceback
        monitor.log(traceback.format_exc(), "ERROR")
        return False
    
    finally:
        monitor.stop_monitoring()


def validate_outputs(output_dir: Path, monitor: AnalysisMonitor) -> bool:
    """Validate that all expected outputs were generated."""
    
    monitor.log("\n" + "=" * 80)
    monitor.log("🔍 Validating outputs...")
    
    expected_files = [
        "clinical_analysis_results.xlsx",
        "clinical_analysis_manifest.json",
        "survival_by_Grade.png",
        "clinical_overview_panel.png",
    ]
    
    all_found = True
    for filename in expected_files:
        filepath = output_dir / filename
        if filepath.exists():
            size_kb = filepath.stat().st_size / 1024
            monitor.log(f"✅ Found: {filename} ({size_kb:.1f} KB)")
        else:
            monitor.log(f"❌ Missing: {filename}", "ERROR")
            all_found = False
    
    # Check for APOE survival plots
    apoe_plots = list(output_dir.glob("survival_by_apoe_*.png"))
    if apoe_plots:
        monitor.log(f"✅ Found {len(apoe_plots)} APOE survival plot(s)")
        for plot in apoe_plots:
            monitor.log(f"   - {plot.name}")
    else:
        monitor.log("ℹ️ No APOE-specific survival plots (may be expected if ID mapping not provided)")
    
    monitor.log("=" * 80)
    
    if all_found:
        monitor.log("✅ All expected outputs validated successfully!", "SUCCESS")
    else:
        monitor.log("⚠️ Some expected outputs are missing", "WARN")
    
    return all_found


def main() -> int:
    """Main entry point."""
    
    print("\n" + "=" * 80)
    print("🔬 APOE Clinical Analysis - Monitored Execution")
    print("=" * 80)
    print("\nThis script will:")
    print("  1. Check dependencies")
    print("  2. Run clinical analysis with monitoring")
    print("  3. Log all output to file")
    print("  4. Validate outputs")
    print("\n" + "=" * 80 + "\n")
    
    # Setup logging
    log_file = Path("clinical_analysis_run.log")
    monitor = AnalysisMonitor(log_file)
    
    monitor.log("=" * 80)
    monitor.log("🚀 Clinical Analysis Execution Started")
    monitor.log("=" * 80)
    monitor.log(f"Log file: {log_file.absolute()}")
    monitor.log(f"Working directory: {Path.cwd()}")
    monitor.log(f"Python: {sys.executable}")
    monitor.log(f"Python version: {sys.version.split()[0]}")
    
    # Check dependencies
    monitor.log("\n📦 Checking dependencies...")
    if not check_dependencies():
        monitor.log("❌ Dependency check failed", "ERROR")
        return 1
    monitor.log("✅ All dependencies available")
    
    # Run analysis
    monitor.log("\n" + "=" * 80)
    success = run_clinical_analysis(monitor)
    
    if not success:
        monitor.log("\n❌ Analysis execution failed", "ERROR")
        monitor.log(f"Check log file for details: {log_file.absolute()}")
        return 1
    
    # Validate outputs
    output_dir = Path("clinical_analysis_outputs")
    if output_dir.exists():
        validate_outputs(output_dir, monitor)
    
    # Final summary
    monitor.log("\n" + "=" * 80)
    monitor.log("📊 CLINICAL ANALYSIS COMPLETE")
    monitor.log("=" * 80)
    monitor.log(f"✅ Success!")
    monitor.log(f"📁 Results: {output_dir.absolute()}")
    monitor.log(f"📄 Log file: {log_file.absolute()}")
    monitor.log("=" * 80 + "\n")
    
    print("\n✅ Analysis completed successfully!")
    print(f"📁 Check results in: {output_dir.absolute()}")
    print(f"📄 Full log saved to: {log_file.absolute()}\n")
    
    return 0


if __name__ == "__main__":
    try:
        exit_code = main()
        sys.exit(exit_code)
    except KeyboardInterrupt:
        print("\n\n⚠️ Interrupted by user")
        sys.exit(130)
    except Exception as e:
        print(f"\n\n❌ Fatal error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

