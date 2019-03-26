# CNV-Content Tool (CNV-ConTool)


**Automatic overlap search between inputed regions or breakpoints and copy number variants databases**

## Before starting

This tool uses the regions or breakpoint positions provided by the user to localize described variants overlapping the same positions. Involved genes are also studied.

The tool was developed to work as a CGI application. The web acessible version of the application will be available soon.

## Dependencies:
+ Web Server with CGI support and configured to handle CGI Programs.
+ python2
+ [biomart](https://pypi.org/project/biomart/)
+ [openpyxl](https://openpyxl.readthedocs.io/en/stable/) 
+ [cgi](https://docs.python.org/2/library/cgi.html), [cgitb](https://docs.python.org/2/library/cgitb.html)
+ other python libraries as sys, collections, itertools and time


## Usage:

To put this tool up and runing just clone the repository, and move the contents to the cgi-bin folder of your webserver.
Point your browser for:

<pre><code> http://[your_web_server]/cgi-bin/CNV-ConTool.py
</code></pre> 

Where [your_web_server] iths the adress of your webserver.

## The output:

The tool retrives two xlsx files with the overlap results, one for breakpoints and another for CNVs.
Examples of output tables are available in the output directory in this repository.

## License:

GPLv2


## Found a bug?

Or maybe just wanto to drop some feedback? Just open an issue on github!
   