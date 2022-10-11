#!/bin/bash

echo "Copying volcanot plot to the current directory"
cp -r /opt/visualization_plugins/volcano_plot/html_data .
echo "Replacing report.tsv file"
rm ./html_data/report.tsv
cp $1 ./html_data/report.tsv