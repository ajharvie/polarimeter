#
"""

"""

import serial

def isOpen(port):
    try:
        erm = serial.Serial(port, timeout=1)
        
        if erm.read():
            return "ready"
            
        else: 
            return "not_ready"
        
        
    except:
        return "not_connected"
        
def ReadDevice(port):
    
    """
    ReadDevice(Port)
    
    Port - String - name of serial port

    """

    ser = serial.Serial(port)
    ser.flushInput()

    
    while True:
        try:
            ser_bytes = ser.readline()
            reading = float(ser_bytes[0:len(ser_bytes)-2].decode("utf-8"))
            return reading
                
        except:
            print("No")
            break
