# 4DMRI Analysis Pipeline

github.com/StanfordCBCL/ArterysDataProcessing/

User Documentation – Melody Dong, 2/8/22

Essential libraries and software:

- SimVascular
- Python 3.6
- VMTK/VTK
- Deformetrica
- pandas
- pydicom
- numpy
- multiprocessing

\*Suggested folder organization:

PREFIX

- MODEL/
  - sv\_model/ (SimVascular project folder)
  - shape\_analysis/
    - cl\_sv/
    - registration\_output/
    - input.txt
    - mag\&lt;min##\&gt;.vtp
    - mag\&lt;max##\&gt;.vtp
- RESULTS/
- PREFIX\_##.vtk
- MR\_headerinfo.txt
- velocity.npy
- dicom\_4dmr.zip

1. **(Optional)** Extract .tar/.zip file downloaded from Arterys which contains all of the images (.dcm dicom files) for the entire scan, which may include the 4DMRI as well as standard MR scans.

  1. Use imagesorting.py script

    1. _python imagesorting.py -h_ \&lt;- to see help options
    2. _python imagesorting.py True 0.1 \&lt;path to .tar/.zip file\&gt;_

To add the path to the .tar/.zip file, you can drag and drop the file into the terminal. If you want to extract multiple files, you can drag and drop multiple .tar files into the terminal.

1. **Convert Arterys 4DMRI scan to readable and manipulated .vtk files.**

Extract .tar/.zip of dicom images from Arterys of the 4DMRI scan to vti images sorted by magnitude. Map magnitude images onto velocity data from numpy array.

  1. Download .tar/.zip file of all DICOM images and velocity numpy array from 4DMRI scan on Arterys.
    1. To download DICOM images, right click on the phonetic ID from Artery list and select &quot;download DICOM&quot;.
    2. To download velocity numpy array from 4DMRI in Artery, drag the 4DMRI scan to the middle viewing window. On the right panel, go to Actions -\&gt; Measure -\&gt; Flow. There will be a _Correction_ button on the most right where you can apply custom or default phase-offset and background correction and Download Velocity Data.
  2. Use arterys\_mapping.py script to extract files from compressed folder and map with velocity.
    1. _Usage: python arterys\_mapping.py [imageFile] [velFile] [saveFileName] [imageExtracted]_

_imageFile: path to compressed folder of dicom files or parent folder where IMAGE/ has already been extracted (by imagesorting.py)_

_velFile: path to velocity numpy array from Arterys_

_saveFileName: file name prefix to save outputs with (e.g PREFIX\_##.vtk)_

_imageExtracted: True or False if compressed folder already extracted_

_python arterys\_mapping.py /arterys\_mapping.py /home/example/dicom\_4dmri.zip /home/example/velocity.npy PREFIX False_

    1. Output: (if extracting DICOM files) IMAGE/ folder with all dicom images sorted into magnitude dicoms by cardiac frame, MR\_headerinfo.txt containing metadata information about venc, heart rate, etc.; PREFIX\_##.vtk saved to parent directory of IMAGE/
1. Create 3D models of the initial timepoint and the max deformed vessel timepoint using SimVascular

  1. Use the vti&#39;s of the magnitude of the 4DMRI scan in time to manually check the timepoint at which maximum deformation of the vessels occurs. This is a visual check and usually signifies the approximate time of peak systole or slightly afterwards due to a phase lag. Create paths, segmentations, and models of the minimum and maximum deformed timepoint in SimVascular
  2. Create a mesh for each model using at least the resolution of the 4DMR scan. On the left panel, right click each mesh and Export to a mesh complete folder.
  3. Copy the mesh-complete.exterior.vtp from each of the exported mesh folder to a new folder (typically a shape\_analysis folder) and rename it as the mag#.vtp with # corresponding to the timepoint from the vti image it was taken from.
1. Use Deformetrica&#39;s shape analysis registration algorithm to interpolate the 3D shape of the PAs in between the initial and max deformed state as well as from the max deformed state back to diastolic state. These two registrations represent systolic upstroke and the gradual systolic downstroke.

  1. Use the postprocessing\_4dmri.py script and an input.txt file (i.e example\_postprocess\_input.txt) to define the parameters and options for this.

    1. In the example\_postprocess\_input.txt file (or you can save a copy of this to a new text file, e.g. PEA.txt which contains the parameters specific to the PEA patient), change the following parameters:

      1. do\_registration = True (set all other &quot;do\_\&lt;\&gt;&quot; options and the centerline\_already\_extracted to False)

        1. If you will be doing Step 6 in conjunction with this step (recommended), set the do\_volumeremeshing to True
      1. data\_path: set to the directory that contains the solid model surface mesh files (i.e mag1.vtp created in Step 4)
      2. Subject\_id: set to a prefixed used for all outputted files from Deformetrica
      3. Dia\_tmpl: set to the .vtp surface mesh model for the initial timepoint
      4. Sys\_tmpl: set to the .vtp surface mesh model for the maximum deformed timepoint
      5. Total\_timepoints: total number of timepoints in the cardiac cycle (i.e total number of vti\_imgs)
    1. To call this script, go to the environment that contains vmtk and type in the following command

      1. _python \&lt;path to postprocessing\_4dmri.py\&gt; \&lt;path to the input.txt file\&gt;_
1. Remesh the output of Deformetrica&#39;s registration output into volume meshes

  1. Use the postprocessing\_4dmri.py script to create the meshes. Creates volume and surface meshes for each geometry outputted by Deformetrica. This may take some time and consume some of your RAM. Change the following in the input.txt file

    1. Change the variables in the input.txt according to what is outlined in Step 5
    2. Set do\_volumeremeshing = True, all other options for analysis to false
  1. Check that these interpolated models match the image. Load the remeshed models and the vti\_img files into Paraview and visually check how well it matches. You can change the opacity and take slices to check.
1. [Create centerlines](#_Instructions_for_extracting) for each remeshed model outputted in step 6.

  1. The outputdir set in the sv\_centerline\_extraction.py script should be the same as the data\_path set in the input.txt file from Step 5/6
1. Resample the 4DMRI vtk files with the remeshed registration outputs at each timepoint to constrain the velocity data to only the MPA.
2. Extract slices from the 4DMRI resampled PA files and compute quantities of interest for each

  1. Additional helper functions and analysis can be added to this part of the code to do more complex analysis

\*For steps 8-9, use the postprocessing\_4dmri.py script and a modified input.txt file

- In the input file, change the following:
  - Do\_registration = False
  - Do\_volumeremeshing = True
  - Do\_resample4dmri = True
  - Do\_qoiextraction = True
  - Centerline\_already\_extracted = True
  - Data\_path: change to path where solid PA models are located
  - Subject\_id: change to desired prefix name for outputted files
  - Mr\_vtk\_dir: change to path where the outputted vtk files from the Step 2 are located
  - Mr\_subject\_id: change to the prefix used for those vtk files from Step 2
- To call this script, go to the environment that contains vmtk and type in the following command:
  - _python \&lt;path to postprocessing\_4dmri.py\&gt; \&lt;path to the input.txt file\&gt;_

##

##

## **Instructions for extracting centerlines using SimVascular:**

1. Import the surface mesh of the geometries outputted from deformetrica (labelled as \&lt;patientname\&gt;\_\_\&lt;timepoint\&gt;\_surfacemesh.vtk). Do this for each timepoint in the scan (i.e all the surface meshes outputted by deformetrica)

  1. In a simvascular project (preferably the one used to create the min and max models), right click on _Models_ from the SV Data Manager panel (left) and select _ **Import Solid Model** _
  2. Select one of the surfacemesh.vtk models and name it _ **def\_\&lt;timepoint\&gt;** _ (e.g _def\_00_ for the H-MA-CL-009\_J0\_\_00\_surfacemesh.vtk model).
  3. When prompted if you would like to extract faces for the model, select _ **Yes** _. Set the _Separation Angle_ to **15**.
  4. The model should show up with a long list of noname faces in the Modeling tab. Identify and label the caps of the model (i.e MPA, LPA, RPA) by selecting the cap on the model. You can do this by hovering over the cap and pressing &#39;_s_&#39; on the keyboard.

    1. ![](RackMultipart20220708-1-xm5ixx_html_8723b7eec8dc4a99.jpg)

1. Once the cap is highlighted and turns yellow, there should be a face in the Face List that is highlighted. That corresponds to the face that you&#39;ve selected on the model. Relabel that face name to the corresponding cap (label these exactly MPA, LPA, or RPA in all caps).
2. When you have finished labelling all of the caps, select those faces in the Modeling Tab, right click and select _Change Type_. Change the type to cap and press OK.
3. Sort the Face List by **Name**. Select all of the noname\_# faces from the list. In the Face Ops on the right, select Combine. Relabel the combined name to wall and change the cap type to wall.
4. Save

1. Run the simvascular centerline extraction python script (sv\_centerline\_extraction.py located in the ArterysDataProcessing repo).

  1. Change the _mdir, outputdir,_ and _total\_timepts._ The _mdir_ is the Models folder in your simvascular project. The outputdir is where you would like the centerlines to be located. Total timepoints is the number of timepoints in the cardiac cycle.
  2. To run the simvascular centerline extraction python script. You can either run from the SimVascular application or from the terminal. (more details located here: [http://simvascular.github.io/docsPythonInterface.html](http://simvascular.github.io/docsPythonInterface.html))

    1. To run from the terminal, locate the simvascular application and enter the following: \&lt;_path to simvascular\&gt; --python -- sv\_centerline\_extraction.py_
    2. To run from the SimVascular GUI, select the Python button on the top toolbar. Go to Text Editor and load in the sv\_centerline\_extraction.py script, and hit run).

Alternative Method for SimVascular Centerline Extraction:

From the terminal SV API:

[directory to simvascular application] --python -- sv\_automated\_cap\_cl.py --path-name=&#39;[path to model&#39;s geometry shape\_analysis/registration\_output]&#39; --prefix-name=&#39;[prefix naming]&#39; --timepoints=[number of timepoints] --clip-distance=0.1 --clip-width-scale=4 --mesh-scale=0.3 --remesh-scale=0.25