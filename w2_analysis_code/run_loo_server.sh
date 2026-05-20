#!/bin/bash

# Script to run LOO_analysis.R persistently on server
# This script ensures the analysis continues running even when your computer turns off

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="$SCRIPT_DIR/LOO_analysis.R"
LOG_FILE="$SCRIPT_DIR/loo_analysis.log"
PID_FILE="$SCRIPT_DIR/loo_analysis.pid"

# Check if R script exists
if [ ! -f "$R_SCRIPT" ]; then
    echo "Error: LOO_analysis.R not found in $SCRIPT_DIR"
    exit 1
fi

# Function to start the analysis
start_analysis() {
    echo "Starting LOO analysis at $(date)"
    echo "Log file: $LOG_FILE"
    echo "Process will continue running even after disconnect..."

    # Start R script with nohup to detach from terminal
    nohup Rscript "$R_SCRIPT" > "$LOG_FILE" 2>&1 &

    # Save process ID
    echo $! > "$PID_FILE"
    echo "Process ID: $(cat $PID_FILE)"
    echo ""
    echo "To monitor progress: tail -f $LOG_FILE"
    echo "To check if running: cat $PID_FILE && ps -p \$(cat $PID_FILE)"
    echo "To stop: kill \$(cat $PID_FILE)"
}

# Function to check if analysis is running
check_status() {
    if [ -f "$PID_FILE" ]; then
        PID=$(cat "$PID_FILE")
        if ps -p "$PID" > /dev/null 2>&1; then
            echo "LOO analysis is running (PID: $PID)"
            echo "Started: $(ps -p $PID -o lstart= 2>/dev/null)"
            return 0
        else
            echo "LOO analysis is not running (stale PID file)"
            rm -f "$PID_FILE"
            return 1
        fi
    else
        echo "LOO analysis is not running"
        return 1
    fi
}

# Function to stop the analysis
stop_analysis() {
    if [ -f "$PID_FILE" ]; then
        PID=$(cat "$PID_FILE")
        if ps -p "$PID" > /dev/null 2>&1; then
            echo "Stopping LOO analysis (PID: $PID)..."
            kill "$PID"
            sleep 2
            if ps -p "$PID" > /dev/null 2>&1; then
                echo "Force killing..."
                kill -9 "$PID"
            fi
            rm -f "$PID_FILE"
            echo "Stopped."
        else
            echo "Process not running, cleaning up PID file"
            rm -f "$PID_FILE"
        fi
    else
        echo "No PID file found, analysis not running"
    fi
}

# Function to show log tail
show_log() {
    if [ -f "$LOG_FILE" ]; then
        echo "Last 20 lines of log:"
        tail -20 "$LOG_FILE"
        echo ""
        echo "To follow live: tail -f $LOG_FILE"
    else
        echo "No log file found"
    fi
}

# Main script logic
case "${1:-start}" in
    start)
        if check_status > /dev/null; then
            echo "LOO analysis is already running!"
            check_status
        else
            start_analysis
        fi
        ;;
    stop)
        stop_analysis
        ;;
    status)
        check_status
        ;;
    restart)
        stop_analysis
        sleep 2
        start_analysis
        ;;
    log)
        show_log
        ;;
    follow)
        if [ -f "$LOG_FILE" ]; then
            tail -f "$LOG_FILE"
        else
            echo "No log file found. Start the analysis first."
        fi
        ;;
    *)
        echo "Usage: $0 {start|stop|status|restart|log|follow}"
        echo ""
        echo "Commands:"
        echo "  start   - Start LOO analysis in background"
        echo "  stop    - Stop running analysis"
        echo "  status  - Check if analysis is running"
        echo "  restart - Stop and start analysis"
        echo "  log     - Show last 20 lines of log"
        echo "  follow  - Follow log in real time"
        echo ""
        echo "Files:"
        echo "  Script: $R_SCRIPT"
        echo "  Log:    $LOG_FILE"
        echo "  PID:    $PID_FILE"
        exit 1
        ;;
esac