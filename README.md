

Citation
--------

If this software is useful for your work, please include one of the following references in your publication or redistribution:

Georg Basler and Zoran Nikoloski (2011): JMassBalance: mass-balanced randomization and analysis of metabolic networks. Bioinformatics, 27(19):2761-2762. https://doi.org/10.1093/bioinformatics/btr448

Georg Basler, Oliver Ebenhöh, Joachim Selbig, Zoran Nikoloski (2011): Mass-balanced randomization of metabolic networks. Bioinformatics, 27(10):1397-1403. https://doi.org/10.1093/bioinformatics/btr145

JMassBalance Installation
--------------------------

Extract the files from the download package to a local directory. Make sure that
the Java Runtime Environment 1.6 or higher is installed by opening a command
line window and typing java -version. The output should look similar to java
version "1.6.0 11". Otherwise, install Java (see http://www.oracle.com).


SBML Support
-------------

For randomizing SBML models, libSBML (http://sbml.org/Software/libSBML) is
required. For Windows, the files libsbmlj.jar and sbmlj.dll are included in the
download package and provide SBML support without the need to install libSBML.

If you wish to use SBML models on Linux or Mac, libSBML has to be in-
stalled with Java support. Note that JMassBalance can also run without SBML
support, so there is no need to install anything if you do not need SBML.
Refer to http://sbml.org/Software/libSBML/docs/cpp-api/libsbml-installation.html
for installing libSBML with Java support on your system.

After the installation, replace the file libsbmlj.jar in the directory where you
extracted the JMassBalance package by the file from the libSBML installation.
On Linux, the latter is usually located at /usr/share/java/libsbmlj.jar.
Additional steps may be required depending on your operating system. Refer
to http://sbml.org/Software/libSBML/docs/cpp-api/libsbml-accessing.html for
more information.


Running JMassBalance
---------------------

JMassBalance can be used via a graphical user interface (GUI) or from the
command line. The GUI is started by opening the file massbalance.jar with the
Java Runtime Environment (java executable). Depending on your system
configuration, this may be done by clicking the file or by executing 

java -jar massbalance.jar

from the command line. For accessing all features of JMassBalance, such as
parallelization, you should execute the program from the command line by typing

java -cp massbalance.jar massbalance/Randomize

in the directory where you extracted the download package. Additional param-
eters are added to specify the file containing the network to be randomized, the
output directory, and other parameters for randomization.

Example:

The download package includes a toy network of the TCA cycle (file 'tca') and
a compounds file with sum formulas (file 'tca.compounds'). Fist, to generate 1000
mass-balanced randomized networks, execute the following command:

java -cp massbalance.jar massbalance/Randomize tca outdir massbalance

The directory 'outdir' will be created, if it does not exist, and the log-files,
mass equivalence classes and randomized networks will be placed therein.

After having generated the 1000 randomized networks, the characteristic path
lengths of the TCA cycle and the randomized networks can be calculated by executing
the following comand:

java -cp massbalance.jar massbalance/Properties tca outdir massbalance pathLength

The results are written to the file 'tca.0-999.massbalance.pathLength' in the
specified directory 'outdir'.

Note that, in order to calculate graph properties, you must specify an index range
which is included in the index range of the previously randomized networks. 


Updates and more information
-----------------------------

For more information, please refer to the Reference Manual
(JMassBalance_manual.pdf). The latest version of JMassBalance can be found at 
http://mathbiol.mpimp-golm.mpg.de/massbalance/.

For accessing the tool from your code, an API Documentation can be found online
and is included in the doc directory of the download package.

For comments and suggestions, please contact Georg Basler.

2023-06-27
