# R-programming-Statistic_analysis

Welcome to my R Assignments Repository! This repository contains a collection of R scripts that I have developed as part of my coursework and projects. Each file in this repository represents a specific assignment or task that I have completed, showcasing my proficiency in R programming and data analysis.
# Channel 4 Analysis Script 

## Overview
This ImageJ macro script, courtesy of Vaishali Jain and Ralph Loring, is designed to process multi-channel TIFF images, specifically CH1, CH2, and CH3 channels. 

  1.This script enables the user to set the scale of the images in pixels and measure the area, mean, standard deviation, etc.

2. Input the ROIs manually for background and cell areas.

3. For each channel (Red, Green, and Blue), open the corresponding TIFF files, unstack them into single color channels, and measure the ROIs.

4. Measure the selected ROIs and label them on the image for visual confirmation.

5. Save modified images with the ROIs and their measurements as TIFF files.

6. Copy and paste the measurement results into an Excel spreadsheet.

7. Clear results and prepare for the next image.

## Why Use This Script?

This script automates multi-channel imaging data analysis, which is an important step in cell-based experiments. It puts together the analysis of images from various color channels that usually represent different markers or sample components. It processes and measures ROIs for each channel to allow accurate counting of cells, area, intensity, and other main features in biological research. This script is perfect for image segmentation and multi-marker analysis, where one wants to segregate a certain signal from different channels, for example, fluorescent markers. This script automatizes ROI selection, measurement, and saving of the results, hence increasing efficiency, reducing human error, and saving time. Besides this, the output will come out in the standard format, easily exportable to Excel for further statistical analysis and reporting, it also stands to be very important in big imaging studies that demand coherence and accuracy in data collection.

## Instructions

Follow these steps to use the script:

1. **Open the Channel 4 (CH4) image in ImageJ**, then click "OK" to start the process.
2. **Enter 5 background ROIs** followed by 't', then **enter 5 cell ROIs** followed by 't', and click "OK".
3. The script will ask you to **copy the results from the "Results" window** to Excel after each channel's result is processed.
4. The script will **automatically process** the red, green, and blue channels and save labeled ROI results for each.
5. Ensure that all channels (CH1, CH2, CH3) are **present in the image directory** as TIFF files before running the script.

## File Structure

- `Channel 4 Analysis Script.ijm`: The ImageJ macro script file.

## Requirements

- ImageJ (or FIJI) installed.
- TIFF files for the red (CH1), green (CH2), blue (CH3), bright field (CH4), and Overlay channels.

## How the Script Works

1. **Set Measurements**: Configures the measurement parameters for the analysis.
2. **Set Scale**: The scale is set to pixels.
3. **Input ROIs**: The user manually selects 5 background and 5 cell ROIs.
4. **Process Channels**: The script processes the CH1, CH2, and CH3 TIFF files (Red, Green, and Blue).
- Each channel is unstacked into its individual color channel.
- The ROIs are measured for each channel.
- The results are copied to the "Results" window and can be saved to Excel.
- Modified images are saved with labeled ROIs.
5. **Clear Results**: Clears the results and prepares for the next image.

