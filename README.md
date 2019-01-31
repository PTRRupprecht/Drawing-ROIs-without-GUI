# Drawing ROIs without GUI
Drawing ROIs with a non-graphical user interface in Matlab

Drawing ROIs is a repetitive task when analyzing calcium imaging data. Despite automatization attempts (e.g., http://neurofinder.codeneuro.org/), most of the time drawing the ROIs manually is still the best solution. This script tries to simplify and streamline manual drawing of ROIs.

### Quick start

To try out the script provided here,

1. Clone the repository locally
2. Open **`Demo_Analysis.m`** in Matlab
3. Follow the notebook instructions and load the (small) demo calcium movie excerpt
4. Once the graphical user interface appears: a) hit the `f` key and draw ROI b) use mouse to navigate.

### The idea: using callback functions without GUI

In order to make this task as easy and intuitive as possible, I used a non-graphical user interface that relies on the user's knowledge of which keys and which mouse buttons to press (inspired by the Google Maps interface). I described the idea of this non-graphical user interface on [my blog](https://ptrrupprecht.wordpress.com/2015/06/24/a-simple-non-graphical-user-interface-in-matlab-keyboard-callback-functions/). The callback functions are defined in [this subscript](https://github.com/PTRRupprecht/Drawing-ROIs-without-GUI/blob/master/non-GUI%20ROI%20analysis/switchImage.m).

[![Drawing ROIs with a non-graphical user interface](nonGIU_ROI_drawing.png)](https://youtu.be/rGTTGCEGvYQ "Drawing ROIs with a non-graphical user interface")

### Overview of key-strokes

`q` : switch to anatomy only\
`w` : switch to dF/F map (response)\
`e` : switch to dF/F map (maximum)\
`r` : switch to anamy overlayed with selected ROIs\
`t` : switch to map of ROIs only\
`y` : switch to map of local correlations

Mouse-click right : full view zoom\
Mouse-click middle : pan (move mouse)\
Scroll-wheel : zoom in and out

`f` : mouse-click left draws ROI outline\
`d` : mouse-click left deletes ROI\
`g` : mouse-click left shifts existing nearby ROI\
`z` : mouse-click left automatically selects a circular ROI\
`x` : mouse-click left automatically selects a ROI based on local correlations\
`1` to `9` : selects ROI size for semi-automatic ROI selection (`x`, `z`)\
`v` : mouse-click left plots a map of local correlations of all image pixels to the selected location (cf. [Junek et al.](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2711456/))

Close figure window : save selected ROIs and timetraces to the Matlab workspace.


### 'Demo_analysis.m' explained

The main script (`Demo_Analysis.m`) loads a 3D stack into the RAM and performs a couple of pre-processing steps before opening the non-graphical user interface (which is basically a simple image opened in Matlab, glued to a couple of key-stroke- and mouse-gesture-specific callback functions). The video linked above gives an overview over the keys (=shortcuts) and how to use them on a real calcium imaging dataset.

To adapt this script to your data, you want to change 'Demo_Analysis.m', not the rest of the code. This program can be used, **as long as your data fit into the RAM**. For larger single calcium recordings, other methods (that use only average images of anatomy or binned videos) must be used.

Variations of this script have been used to extract activity of neurons from calcium imaging data in the following papers:  [paper 1]( http://dx.doi.org/10.1016/j.cub.2017.11.007), [paper 2](https://www.osapublishing.org/boe/abstract.cfm?uri=boe-7-5-1656), [paper 3](https://doi.org/10.1016/j.neuron.2018.09.013).
