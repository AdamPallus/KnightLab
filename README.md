# Knight Lab Analysis Code
Current Dashboard: https://gaze.shinyapps.io/dashboardtest2/
## To update a dashboard with new data:
* Use combineMAT.mat file in Matlab to combine .mat data files from experiment computer and EyeSeeCam
* Put resulting CSV files in a folder called "data" inside this repo
* Run measure_script.R to measure trials and write data file to dashboards/output/ 
* Update the first cell of dashboardtest2.Rmd to read the new file
* Run shiny app locally
* Click publish to publish new dashboard to shinyapps.io
