Individual files can be viewed in file-tree or entire folders downloaded as zip files:  
- `/mat` contains ~15 GB of data  
- `/classes` & `/scripts` contain small .m files only


----------


#### **/classes**  

Matlab object classes with core functions for analysis of tiff files: 

- ##### **SlidebookObj** 
    Superclass, one object per tiff file.  
    Provides base functions for `avp4DmaxObj` - not used directly.  
    Additional documentation: [github.com/bjhardcastle][1] 

- ##### **avp4DmaxObj** 
    Subclass of `SlidebookObj` for polarization experiments.  
    One object per tiff file (single imaging plane or max intensity projection).

- ##### **avp4DsuperObj** 
    Independent class for groups of associated `avp4DmaxObj`  
    One object per recording.  
    Contains multiple objects:  
    - active channel (GCaMP) 
      - MIP (single `avp4DmaxObj`)
      - individual layers from stack (`[1xN]` `avp4DmaxObj` array)
    - static channel (tdTomato) 
      - MIP (single `avp4DmaxObj`)
      - individual layers from stack (`[1xN]` `avp4DmaxObj` array)

- ##### **/additionalfuncs, /icons** 
    Extra functions for plotting or viewing data, many third-party, including modified version of James Strother's Neuron Image Analysis GUI for playback of tiff files: [bitbucket.com/jastrother][2] 

----------


#### **/mat**

Matlab .mat files containing data and object arrays: 
- ##### **/aux_data_backup**  (~5 GB)  
    Minimal storage of data.  
    Subfolders organized by neuropil, cell-type, Gal4 line, then date.  
    Files associated with a single recording are grouped, including max intensity projections, ROIs, masks, etc. 

- ##### **/objarrays**  (~10 GB)
    One .mat file per cell-type.  
    Each contains a `[1xN]` `avp4DsuperObj` array  

- ##### **/plotting, /psi**  (<20 MB)
    Additional derived data stored to save processing time.

----------


#### **/plots**
Printed plots from Matlab: .pdf & .png organized by figure.  
 - `commandWindow.log` contains record of stats returned during analysis when run via `make_all_figs.m`

----------


#### **/scripts**
- ##### **/load_object_functions**
    Loads .mat files in objarrays along with some additional variables/processing.

- ##### **/make_figure_panels** 
    Scripts for running analysis and printing plots.

- ##### **/plotting** 
    General functions for visualizing data or specific scripts to generate one-off plots.  
- ##### **/refresh** 
    Scripts to re-analyze and store .mat files in `/mat/plotting/` 

----------

#### **explore data**
All folders are required, with the structure preserved.

Add certain folders to Matlab's search path by running `pathsAVP.m` from the root folder.  

Load an object array with one of the files in `/load_object_functions`: 

    loadDmDRA   % loads [1xN] avp4DsuperObj array as variable 'x'
  
A quick way to explore selectivity and tuning values in each recording is through a GUI for making masks/ROIs:

    maskLayers(x(1))    % open GUI with the first recording in the array

@[osf](x7asg)

Selectivity and tuning maps can be toggled on/off. Scatter plots can also be produced to examine polarotopy within the mask (applies threshold in 'selectivity' text box):
@[osf](vgzyu)

  [1]: https://github.com/bjhardcastle/SlidebookObj
  [2]: https://bitbucket.org/jastrother/neuron_image_analysis/src
