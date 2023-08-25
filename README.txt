~~~ CS661A ASSIGNMENT 4: SAMPLING AND RECONSTRUCTION ~~~

~ Installation (optional) ~
I have included a requirements.txt file for installing relevant libraries. Skip the installation if all relevant libraries are installed already on your system. See instructions here on how to use requirements.txt: https://stackoverflow.com/questions/7225900/how-can-i-install-packages-using-pip-according-to-the-requirements-txt-file-from

~ Running the program ~
Make sure sampling.py and Isabel_3D.vti are in the same directory. Then open the terminal in the same directory and enter the following:

$ python3 sampling.py --percent {sampling percentage} --method {reconstruction method}

Here is an example. Let's say we want to do 3% sampling and perform reconstruction using the linear method. Then the command will be:

$ python3 sampling.py --percent 3 --method linear

Change the command above to run the program with your desired sampling percentage and reconstruction method.

~ Seeing the program outputs ~
After the program has run, two files will be written. By default, the sampled points will be saved in a file named sample.vtp and the reconstructed data will be saved in a file named {reconstruction method}_recon.vti. For example, if you choose nearest interpolation method, then the file will be named nearest_recon.vti. If you choose linear interpolation method, then the file will be named linear_recon.vti. Note that running the program again will overwrite the previous files.
