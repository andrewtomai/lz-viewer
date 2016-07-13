LocusZoom Association Viewer
============================

Summary
=======
The LocusZoom Association Viewer allows the user to easily view the results of their own genome wide association tests, including PLINK, RAREMETALWORKER, and EPACTS.  Upon runnning this program, the user browses to localhost (127.0.0.1) on their browser, at the port they specified. Three different types of graphs are available:  

* Manhattan Plot
* LocusZoom Plot
* QQ Plot
	* Currently, only available for EPACTS files
 
Command Line Options
--------------------
The only required command line argument is the filename, which should be directly after the program name.
The program is ran with multiple options including:

--help         list off all of the command line options
--port         you can specify a port at which to view the results
--range        specify a specific range of format [CHROMOSOME]:[START POSITION]-[END POSITION] to view.

	* if range is specified, the default homepage is the LocusZoom Plot with specified range.
	* if range is not specified, the default homepage is the Manhattan Plot

Sample use of command line options:
-----------------------------------
lz_assoc_viewer my_results.gz 11:113580000-114150000 --port 5541 


Manhattan Plot
--------------
The manhattan plot is a scatter plot that plots all points.  However, if the point is not significant enough, it is binned.  Points that are significant enough have tooltips and are clickable.  Upon click, the user is taken to a LocusZoom Plot with a region that is centered on the point clicked.

QQ Plot
-------
Eventually, the QQ plot will include clickable tooltips just like the manhattan plot.

LocusZoom Plot
--------------
The LocusZoom Plot is a scatterplot of a specific region of a chromosome.  Points have tooltips that show additional information.

