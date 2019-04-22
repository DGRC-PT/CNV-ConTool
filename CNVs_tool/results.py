#!/usr/bin/env python
#coding=UTF-8

import sys
sys.path.append('/home/analise_liWGSseqs/software/biomart')
sys.path.append("/usr/lib/pymodules/python2.7/openpyxl")
from biomart import BiomartServer
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill
from openpyxl.styles.borders import Border, Side
import time
import cgi
import cgitb
import itertools
from collections import OrderedDict
import cross_deletions_with_dgv_results_V8
cgitb.enable()

def index():
        """ Show the default page
        """
        print 'content-type: text/html'
        print #Blank line to divide headers
        print '</html>'

index()


def indexx():
	print 'Content-type: text/html; charset=utf-8'
	print """
 <html>
 <meta name="viewport" content="width=device-width, initial-scale=1.0">
<link href="https://fonts.googleapis.com/css?family=Arial" rel="stylesheet">
<style>
* {
  box-sizing: border-box;
}

body {
  background-color: #f1f1f1;
}

#regForm {
  background-color: #ffffff;
  margin: 80px auto;
  font-family: Arial;
  padding: 40px;
  width: 70%;
  min-width: 300px;
}

hgh {
  text-align: center;
  font-family: Arial;
}

aaa {
  padding:40px;
  font-family: Arial;
}

h1 {
  text-align: center;
  font-family: Arial;
}

</style>


<head>
<center><title>CNV-ConTool</title></center>
</head>
<aaa>
<center><img src="https://cld.pt/dl/download/c3ba7b53-2d29-46fa-a2e0-c6dfbf2310e9/banner_HMSP.jpeg" width=622 height=153 border=0 alt=""></center>
<center><h1>CNV-Content Tool - Search Results</h1></center>
<center><p>The retrieved data is compiled in two xlsx files, one for breakpoint and another for CNV results, composed by a summary table, and a set of more complete tables, one for each database used.</p></center>
<br>
<br /><hgh><center><b> While you wait for the results, please fill out our <a  href="https://goo.gl/forms/WKiyDJtXbXuSgDYG2" target="_blank">usage survey</a>. Thank you! </center></b></font><br>
<br /><center>Results will appear shortly...</center>
</aaa>
</html>
        
"""

indexx()
cgitb.enable()
form=cgi.FieldStorage()
#print(form)
cross_deletions_with_dgv_results_V8.indexx()
cross_deletions_with_dgv_results_V8.exect(form)
cross_deletions_with_dgv_results_V8.remake()
