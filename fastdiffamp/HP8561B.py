# GPIB data grabber for HP8561B Spectrum Analyzer. This script uses the Prologix GPIB-to-USB
# controller to connect to the HP-IB interface using a serial bus.
# When the script is run Python connects to comport 6 (defined in varible "comport"). The
# GPIB-to-USB controller is the set to communicate with GPIB-address 18 (defined in varible 
# "GPIB_addr"). The data corresponding to what is shown on the HP8561B sceen. If varible "plot" is 
# True, the data is plotted.

import serial
import time
import decimal
import matplotlib.pyplot as plt

# Function that sents a command to the HP8561B and waits for an answer.
def GPIB_ask(cmd_str, bytes):
	GPIB.write(cmd_str)
	GPIB.write("++read eoi\n")
	answer = GPIB.readline()
	return answer

def write_data(x_val, y_val, N):
	filename = "GPIB_data_" + time.strftime("%y%m%d%H%M") + ".csv"
	GPIB_id = GPIB_ask("ID?;\n",16)
	with open(filename, "w") as myfile:
		myfile.writelines( GPIB_id[0:-1] + " data file.\n")
		myfile.writelines("Frequency [Hz], Amplitude [dBm]\n")
		for i in range(N-1):
			myfile.writelines(str(x_val[i]) + ", " + str(y_val[i]) + "\n")
		myfile.writelines(str(x_val[N-1]) + ", " + str(y_val[N-1]))
	myfile.close()
	print("Data written to", filename)

if __name__ == '__main__':
	comport = 8		# Comport for the GPIB-to-USB controller
	GPIB_addr = 18	# GPIB address of the HPIB.
	plot = True		# If True, the data will be plotted.

	GPIB = serial.Serial( comport-1, 9600, timeout=10)	# Initialise the serial port.

	GPIB.write("++mode 1\n")						# Set controller mode.
	GPIB.write("++addr " + str(GPIB_addr) + "\n")	# Set GPIB address.
	GPIB.write("++auto 0\n") 						# Disable read-after-write.

	GPIB_id = GPIB_ask("ID?;\n", 16)
	print("Getting data from ", GPIB_id)

	GPIB.write("TDF P;\n")				# Formats trace data to ASCII decimal values.
	data = GPIB_ask("TRA?;\n", 10240)	# Asks for the amplitudes of trace A.

	#time.sleep(10)
	amplitudes = data.split(',')
	n = len(amplitudes)
	print("Number of points: ", n, "\n")

	frequency = []
	start_freq = GPIB_ask("FA?;\n",16)	# Asks for the start frequency.
	end_freq = GPIB_ask("FB?;\n",16)	# Asks for the stop frequency.

	for i in range(n):	# Loop that converts amplitude from datatype string to float 
						# and creates an array of frequency values corresponding to 
						# start_freq and slut_freq.
		amplitudes[i] = float(amplitudes[i])
		frequency.append( float(start_freq)+(i)*(float(end_freq)-float(start_freq))/(n-1) )


	#print "Wavelength: ", len(frequency)
	#print "Amplitudes: ", len(amplitudes)

	write_data(frequency, amplitudes, n) 	# Saves data.
	print("Done.")

	GPIB.write("++loc\n")	# Return HP70950 to local control.
	GPIB.close();			# Closes the serial port.

	# Plots the data.
	if plot:
		plt.plot(frequency, amplitudes)
		plt.xlabel('Frequency [Hz]')
		plt.ylabel('Amplitude [dBm]')
		plt.show()
