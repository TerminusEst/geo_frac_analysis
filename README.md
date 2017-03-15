#**Geological Fracture Analysis**
=================================

A Python module for the spatial analysis of geological fractures from shapefile data. Includes plotting functions.

This module takes shapefiles containing either polylines or points and spatially analyses them by user defined inputs
(cell size in km and angle bins).

Analysis outputs include:

1. Number density.

2. Length density.

3. Number anisotropy.

4. Group Dominance Frequency (GDF).

Single or multiple analysed shapefiles can then be put into a fancy plot with rose diagrams like so:
![First](https://cloud.githubusercontent.com/assets/20742138/17325319/d67d1608-58a2-11e6-8e6b-59bd155cb080.png) "Ireland, 100k Fracture Map"

## **Quick Code Example**
Import the module:

```python
>>> import geo_frac_analysis as GFA
```

Set address for shapefile which contains polyline fractures.

```python
>>> address1 = "/Data/wexfaults.dbf"
```

Set the cell size to be 10km, and angle bins per 180 degrees to be 3 (i.e., 30 degree bins):

```python
>>> cell_size, angle_bins = 10, 3
```

Create the FracAnalysisPoly object:

```python
>>> a = GFA.FracAnalysisPoly(address1, cell_size, angle_bins)
```

Plot the analysed data with number density patches:

```python
>>> GFA.FancyPlot(a, patches = "Number")
```

##**Installation**
To install using pip:

```python
pip install geo_frac_analysis
```

##**Dependencies**
pyshp

matplotlib

numpy

##**Author**
Written by Sean Blake in Trinity College Dublin, 2014-2016

Email: blakese__at__tcd.ie

GITHUB: https://github.com/TerminusEst

Uses the MIT license.

##**Detailed Use**

Using **a** as our analysed polyline object (see quick example above), the basic plot

`>>> GFA.FancyPlot(a)` gives the following plot:
![2nd](https://cloud.githubusercontent.com/assets/20742138/17325972/6819e1f6-58a6-11e6-8a93-69a8fd4405fb.png) "Basic Plot"

To include squares (size = cell_size) showing number density, we can call:

`>>> GFA.FancyPlot(a, Patches = "Number")` Patches can also equal "Length" or "NumberAnisotropy"
![3rd](https://cloud.githubusercontent.com/assets/20742138/17326097/119ca132-58a7-11e6-9c67-be9c4ac84d89.png) "Basic Plot + Patches"

An example of the plot working with all of the bells and whistles:

`>>> GFA.FancyPlot(a, Rose = True, Fractures = True, Patches = "NumberAnisotropy", Circles = True, SquareNumbers = True, FigureNumber = 1)`
![4th](https://cloud.githubusercontent.com/assets/20742138/17326180/76581b24-58a7-11e6-8206-a046e9338b2b.png) "Fancy Plot"

Suppose we had another shapefile of a different region, and we wanted to plot them together. We can simply do the following:

```
>>> address2 = "/dungarv.dbf`
>>> b = GFA.FracAnalysisPoly(address2, cell_size, angle_bins)
>>> analysed_list = [a, b]
>>> GFA.FancyPlot(analysed_list, Rose = True, Fractures = True, Patches = "NumberAnisotropy", Circles = True, SquareNumbers = True, FigureNumber = 1)
```
![5th](https://cloud.githubusercontent.com/assets/20742138/17328810/9b4708d6-58b7-11e6-9835-451495146d92.png) "Fancy Plot"

Shapefiles which contain points can also be added:
```
>>> address3 = "/mypointdata.dbf"
>>> p = GFA.FracAnalysisPoint(address3, cell_size ,angle_bins)
>>> analysed_list.append(p)
>>> GFA.FancyPlot(analysed_list, Rose = True, Fractures = True, Patches = "Number", Circles = False, SquareNumbers = False, FigureNumber = 1)
```
![6th](https://cloud.githubusercontent.com/assets/20742138/17328979/5db2638e-58b8-11e6-893e-b0870d75e7bd.png) "Fancy Plot"
**NOTE** No length analysis can be undertaken on a group of analysed shapefiles if one or more of them contain points.

If you want to plot one rose diagram per region, we use the second plotting function:

`GFA.FancyPlotTotals(analysed_list, Fractures = True, Circles = True, FigureNumber = 1)`
![7th](https://cloud.githubusercontent.com/assets/20742138/17329097/f150c540-58b8-11e6-9707-2fd6b9bad04c.png) "Fancy Plot"

Now that you are happy with the analysis, you can save the analysed data from **a** into a shapefile. This can be done as follows:

```
>>> output_address = "/my_output"
>>> a.save_output(output_address)
```
Once this is done, you will have a shapefile with fields:
```
[('DeletionFlag', 'C', 1, 0),
 ['X', 'F', 10, 8],
 ['Y', 'F', 10, 8],
 ['N0', 'F', 10, 8],    
 ['N1', 'F', 10, 8],
 ['N2', 'F', 10, 8],
 ['L0', 'F', 10, 8],
 ['L1', 'F', 10, 8],
 ['L2', 'F', 10, 8],
 ['N_tot', 'F', 10, 8],
 ['N_Anisotro', 'F', 10, 8]]
```

where:

X, Y = X and Y coordinates of squares

N0-Nn = number density per bin (bins 0-n)

L0-Ln = length density per bin (bins 0-n)

N_tot = true number density per square (so as to avoid double counting fractures which are in >1 angle bins)

N_Anisotro = number anisotropy per square


