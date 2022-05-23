import subprocess
import sys
import time

if __name__ == "__main__":  # a guard from unintended usage
    proc = subprocess.Popen(["./real-time-lidar"],  # start the process
                            stdin=subprocess.PIPE,  # pipe its STDIN so we can write to it
                            stdout=sys.stdout, # pipe directly to the output_buffer
                            universal_newlines=True)

    time.sleep(0.5)
    print("11", file=proc.stdin, flush=True)
    time.sleep(0.1)
    print("n", file=proc.stdin, flush=True) 
    time.sleep(0.1)
    print("", file=proc.stdin, flush=True)  

    time.sleep(1)
    proc.kill()

