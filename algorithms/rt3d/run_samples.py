import subprocess
import sys
import time
import threading

def thread_function(proc): 
    print('in thread function')   
    time.sleep(2.5)
    print('in thread function2')   
    proc.kill()
    print('in thread function3')   
    

if __name__ == "__main__":  # a guard from unintended usage
    
    for i in range(1,2):
    
        proc = subprocess.Popen(["./real-time-single-photon-lidar/RT3D.exe"],  # start the process
                                stdin=subprocess.PIPE,  # pipe its STDIN so we can write to it
                                stdout=subprocess.PIPE, # pipe directly to the output_buffer
                                universal_newlines=True)

        time.sleep(0.5)
        
        print(i, file=proc.stdin, flush=True) # Choose dataset
        time.sleep(0.1)
        print("", file=proc.stdin, flush=True) # Default SBR
        time.sleep(0.1)
        print("1", file=proc.stdin, flush=True) # Algorithm 1
        time.sleep(0.1)
        print("y", file=proc.stdin, flush=True) # Set hyperparameters manually
        time.sleep(0.1)
        print("1", file=proc.stdin, flush=True) # Scale ratio    
        time.sleep(0.1)
        print("1", file=proc.stdin, flush=True) # Number of points per pixel    
        time.sleep(0.1)
        print("", file=proc.stdin, flush=True) # Upsampling rate
        time.sleep(0.1)
        print("", file=proc.stdin, flush=True) # Intensity threshold    
        time.sleep(0.1)
        print("", file=proc.stdin, flush=True) # Proportion of dilated intensity
        time.sleep(0.1)
        print("", file=proc.stdin, flush=True) # Number of iterations
        time.sleep(0.1)
        print("2", file=proc.stdin, flush=True) # Neighbourhood size
        time.sleep(0.1)
        print("", file=proc.stdin, flush=True) # APSS kernel size
        time.sleep(0.1)
        print("", file=proc.stdin, flush=True) # Regularization par for background
        time.sleep(0.1)
        print("", file=proc.stdin, flush=True) # Max reflectivity
        time.sleep(0.1)
        print("", file=proc.stdin, flush=True) # Gradient reflectivity step size
        time.sleep(0.1)
        print("", file=proc.stdin, flush=True) # Depth step size
        time.sleep(0.1)
        print("", file=proc.stdin, flush=True) # Upsampling rate
        
        t = threading.Thread(target=thread_function, args=(proc,))
        t.start()
        
        while True:
            data = proc.stdout.readline()   # Alternatively proc.stdout.read(1024)
            print(len(data))
            sys.stdout.write(data) 
            if len(data)==0:
                break

