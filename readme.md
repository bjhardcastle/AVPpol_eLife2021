Two-photon calcium imaging  analysis code for investigation of polarization-sensitive neurons in the Drosophila anterior visual pathway (AVP)

Corresponding data can be found in the [project repository][3] at the Open Science Framework.

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

#### **notes**

Convert polarizer angles as-recorded to convetion used in publication: `wrapTo180(-theta-270)`

---------

  [1]: https://github.com/bjhardcastle/SlidebookObj
  [2]: https://bitbucket.org/jastrother/neuron_image_analysis/src
  [3]: https://osf.io/3tsd6/
 
