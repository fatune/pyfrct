pyfrct
======

*pyfrct* is a python tool for basic manipylation on a structural geological data
(attitudes of a plane, dir angel, dip angle etc) and for plotting standard 
stereonet diagram.

Using pyfrct
------------

To use pyfrct to plot standard stereonet diagramm with density gray-shaded map 
of the fractions on the backround simply tipe following:

```
# compute density grid of a fracture planes
python density.py -s 5 -r 100 -i data.txt -o data.density

# plot diagram and store it to data.pdf
cat data.txt | python fr_stereo.py -d data.density -o data.pdf -u -r
```


