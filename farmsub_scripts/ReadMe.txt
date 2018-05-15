This is a How-to for the submission of g4sbs jobs *on the JLab batch farm*.

Pretend you want to submit 1M DIS events in the nDVCS configuration.
First of all, you need to build a .mac script to define your configuration (exp. type, spectrometer angles, etc...), and your generator (type, min/max polar azimuthal angle, number of events, etc). 
For the case we consider in this example, the .mac script is already written: ndvcs_dis.mac
To modify something in it, or create a new one, refer to the g4sbs documentation at the following link:

https://hallaweb.jlab.org/wiki/index.php/Documentation_of_g4sbs

Once you have your .mac script, you can call the corrsponding python script (which in this case is makejobs_dis.py), with the following command:

> python2 makejobs_dis.py ndvcs_dis.mac 50

(in our case, the .mac script generates 20000 DIS electrons, one needs to throw 50 jobs to get 1M evts total).

If for instance you want to turn off the scintillation process (with preinit script preinit_ckov_noscint.mac), add the name of the script at the end:

> python2 makejobs_dis.py ndvcs_dis.mac 10 preinit_ckov_noscint.mac 

There is, on top of each .py file, a few lines which recall this.

Note: *it is very important that the output file name defined in the .mac file (with command /g4sbs/filename) is the same that in the output file name called at line 52 of the .py file* otherwise your output file will be lost and your computation time will go down the drain.

Note 2: For the sake of being a good batch farm citizen, I would always advise, before running large stats, to throw 1 test job (low stats, etc...) to check that you have some output at the end, and that this output corresponds to what you expect (i.e. that you have the variables you expect, etc...). 
