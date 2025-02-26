A first rough explanation:

This repository provide a tool-box (in **python** and **C++**) for composite laminate FE mesh creation, mainly using [gmsh](https://gmsh.info/) (for mesh creation) and [Eigen](https://Eigen.tuxfamily.org/index.php?title=Main_Page) (for **C++** Matrix classes) libraries.\
The VTK format is mainly used so the user could use [paraview](https://www.paraview.org/) for visualisation.\
The objective is to provide inputs for [Abaqus](https://www.technia.co.uk/software/simulia/abaqus/) and [Dune-composites](https://gitlab.dune-project.org/anne.reinarz/dune-composites).

* **Csection_mesh** contains **python** scripts which generate the input file (*.geo*) for [gmsh](https://gmsh.info/).

* **Csection_mesh** for Csection part but also flat laminate (depending on the option chosen).

* **C++** folder contains tools for FE mesh manipulation (grid transformation, local material orientation attribution, etc).

* All **C++** tools use the same headers libraries locate in **C++/shared_headers**.
* **C++/gridModification** is used for material orientation attribution and potentially grid transformation (a ramp for the Csection or a wrinkle).\
* **C++/odb2vtk** provides examples of ODB results reader (**C++/odb2vtk/Abaqus2txt_examples**) to extract Abaqus results as a list save in a *.txt* file; and a **C++** function (**C++/odb2vtk**) used to associate these extracted data to the corresponding mesh in VTK format.


* Ini files are read with the [mINI](https://github.com/metayeti/mINI) library within C++ codes.

---

Procedure to compile any C++ code in this repository:

For windows:
all cmake command have to be followed by:
```
-G "Unix Makefiles"
```
You might need to install make and cmake and add them to path

---

* [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) installation (as provided by Eigen doc) and linkage to these tools:

```
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
git checkout 3.3
mkdir build_dir
cd build_dir
cmake ..
```
You may need to specify your compiler.
```
cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc
make install
```
This last needs administrator priviledge.\
Now you have to specify the path_to_eigen in your environment variable at the name *EIGEN_INC*.

* **gridModification** and **odb2vtk**:
```
cd build
cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_C_COMPILER=gcc
make
```

---\
---
* **APPLICATION**

For now the python folders (**Csection_mesh** and **Stiffener_mesh**) are used as working directory so a commun use of the tools is:
```
cd Csection_mesh
python write_parameters.py
```
You can now check *parameters.txt* and modify it if needed.
```
python main.py
```
*main.py* provide the command line to use [gmsh](https://gmsh.info/).

Output: 
* A mesh in .msh2 format.

From now, you may need to mkdir an **Abaqus** and/or **DUNE** folders if they are not in the current folder. Then:
```
../C++/gridModification/build/gridMod
```
Output:
* A result in .vtk format.
* .inp files in **Abaqus/** for Abaqus
* .msh and *_ori.txt files in **DUNE/** for Dune-composites




---
**Precedure to test the delam in ramp/corners**

``` cd /path/to/Csection_mesh```
``` python write_parameters_for_test_delam_function.py```
* in write_parameters_for_test_delam_function.py:
    * See l28 for geometry changes
    * l29: the interlayer where the delam should be (in real part) is selected default: 20
    * l96: The discretization of the Z dir is set to only disretize the central zone (where the ramp is)


``` python main.py```



```(cd ../C++/gridModification/build && make)```
* in gridTransformation.h:
    * The functions ```TESTtoFlatCoordinateSys``` (l1240) and ```TESTinverse_ramp_gridTransformation``` (l1427) replicate the ones available in DUNE composites DefectParamaters.hh
* in gridMod.h:
    * The function ```TESTinitialiseDamage``` (l1080) does the loop over elements to test if the element belong or not to the delam area
    * At the end of this function (l1148) an extra loop is operated to transform (flatten) the mesh, if the booleen ```config['GENERAL']['flatten(b)']``` is true in the ini file.


```../C++/gridModification/build/gridMod```


* The ```Cspar_test_delam.vtk``` will be written.
    * In paraview:
        * Do threshold on scalar:**marker**, selected component:**X** with the value ```lower threshold = upper threshold = 2``` (corresponding to the interlayer)
        * Plot the variable **delam**
