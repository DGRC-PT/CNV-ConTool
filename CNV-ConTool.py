#!/usr/bin/env python
#coding=UTF-8

import cgi
import os
import cgitb
from time import time, localtime, strftime
import datetime
import calendar

cgitb.enable()

clock=strftime("%a, %b %d %Y %H:%M:%S", localtime())

def index():
        """ Show the default page
        """
        print 'content-type: text/html'
        print #Blank line to divide headers
        print '</html>'

index()



def showForm1():
	"""Show a form
	"""
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

aaa {
  font-family:Arial;
  font-size:20px;
  line-height: 0.7;
}

tyt {
  font-family:Arial;
  font-size:25px;
}

rodape {
  font-family:Arial;
  font-size:17px;
}

#regForm {
  background-color: #ffffff;
  margin: 20px auto;
  font-family: Arial;
  font-size:20px;
  padding: 20px;
  width: 70%;
  min-width: 100px;
  line-height: 0.7;
}

h1 {
  text-align: center;  
}

input {
  padding: 10px; 
  font-size: 20px;
  font-family: Arial;
  border: 1px solid #aaaaaa;
  line-height: 0.5;
}


head2 {
  padding: 10px;
  width: 100%;
  font-size: 17px;
  font-family: Arial;
  border: 1px solid #aaaaaa;
}

/* Mark input boxes that gets an error on validation: */
input.invalid {
  background-color: #ffdddd;
}

/* Hide all steps by default: */
.tab {
  display: block
}

.nn {
  display: none;
}

button {
  font-family: Arial;
}

button:hover {
  opacity: 0.8;
}

.column {
  width: 200px;
  margin: 5px;
  display: inline-block;
  text-align: left;
}
.container {
  margin:auto;
  text-align: center;
}

#prevBtn {
  background-color: #bbbbbb;
}

input[type=checkbox]
{
  /* Double-sized Checkboxes */
  -ms-transform: scale(1.2); /* IE */
  -moz-transform: scale(1.2); /* FF */
  -webkit-transform: scale(1.2); /* Safari and Chrome */
  -o-transform: scale(2); /* Opera */
  padding: 10px;
}

/* Make circles that indicate the steps of the form: */
.step {
  height: 15px;
  width: 15px;
  margin: 0 2px;
  background-color: #bbbbbb;
  border: none;  
  border-radius: 50%;
  display: inline-block;
  opacity: 0.5;
}

.step.active {
  opacity: 1;
}

/* Mark the steps that are finished and valid: */
.step.finish {
  background-color: #4CAF50;
}
</style>
<br><br />
<center><img src="https://cld.pt/dl/download/c3ba7b53-2d29-46fa-a2e0-c6dfbf2310e9/banner_HMSP.jpeg" width=622 height=153 border=0 alt=""></center>
<head>
<center><title>CNV-ConTool</title></center>
</head>
<head2>
<aaa>
<center><h1>CNV-Content Tool (CNV-ConTool)</h1></center>
<center><p><b>Automatic overlap search between inputed regions or breakpoints and copy number variants databases</b></p></center>
<center><p>This tool uses the regions or breakpoint positions provided by the user to localize described variants overlapping the same positions. Involved genes are also studied</p></center>
<center><p> The tool retrives two xlsx files with the overlap results, one for breakpoints and another for CNVs.</p></center><br>
</head2></aaa>


<body>
<form id="regForm" method=POST action="CNVs_tool/results.py">
  <!-- One "tab" for each step in the form: -->
  <div align="left">
  <p><b><tyt>&emsp;&emsp;General options</tyt></b></p>
  </div> 
  <hr>
  <div class="container"><b> <p>Reference genome assembly:</b></p>
	<div class="column"><center><input type="radio" value="hg19" name="version">Hg19 </div>
	<div class="column"><p><input type="radio" value="hg38" name="version" checked="checked">Hg38 </p></center></div>
  </div>
  <div class="container"><b> <p>Databases:</b></p>
	<div class="column"><p><center><font color="grey" >Type of alteration:</font></p>
		<p><input type="radio" name="tt" value="loss">Loss</p>
		<p><input type="radio" name="tt" value="gain">Gain</p>
		<p><input type="radio" name="tt" value="inversion">Inversion</p>
		<p><input type="radio" name="tt" value="all" checked="checked" >All</p><br><br></center>
	</div>
	<div class="column"><p><font color="grey" >Databases:</font></p>
		<p><input type="checkbox" name="dats[]" value="DGV" >DGV</p>
		<p><input type="checkbox" name="dats[]" value="1000Genomes" >1000Genomes</p>
		<p><input type="checkbox" name="dats[]" value="ClinGen" >ClinGen</p>
		<p><input type="checkbox" name="dats[]" value="ClinVar" >ClinVar</p>
		<p><input type="checkbox" name="dats[]" value="CoeCoop" >Coe & Cooper</p>
  </div>
  <div class="container"><b> <p>Overlap cutoff (1-100)%:</b></p>
    <p><input type="text" name="perc" value="70"></p></center><br>
  </div>
  <div class="tab"><center><b>Genes flanking region size:</b>
  <p><font color="grey" > Added to each side of the gene </font></p>
  <p><font color="grey" > The region size can be given as a percentage (calculated relativelly to gene size, followed by %)</font></p>
  <p><font color="grey" > or an absolute value (followed by the unit of mesurement bp and kb are accepted). </font></p>
    <p><input type="text" name="percgene" placeholder="ex.: 70% or 100bp or 200kb"></p></center>
  <div class="tab" ><b>Select type of alteration to be analysed:</b>
    <p><select name="ttt" form="regForm" onchange="yesnoCheck(this);"></p>
      <option value="">Selection options</option>
      <option value="breakpoint">Breakpoint</option>
      <option value="cnv">CNVs</option>
      <option value="bpcnv">Breakpoint & CNVs</option>
      </select>
  </div>
      <div id="ifmeh" style="display: none;" align="left"> <p><b><tyt>&emsp;&emsp;Breakpoint analysis</tyt></b></p>
		<hr>
		<center><b>Input breakpoint positions:</b>
		<p><font color="grey" > One per line. Examples of formats with eligible formats are given in the text area bellow</font></p>
		<p><textarea style="font-size: 17px" name="message" rows="5" cols="50">
chr16:174000-174500
t(16;17)_chr16:174000-174500 (advised format)</textarea></p>
		<b><p>Breakpoint flanking region:</b></p>
		<p><font color="grey" > Region added to each side of the breakpoint (in bp). </font></p>
		<p><input type="text" name="bpflank" value="100"></p></center><br>
      </div>
      <div id="ifins" style="display: none;" align="left"> <p><b><tyt>&emsp;&emsp;CNV or specific genomic region analysis  </tyt></b></p>
		<hr>
		<center><b>Input genomic regions or CNVs:</b>
		<p><font color="grey" > One per line. Examples of eligible formats are given in the text area bellow</font></p>
		<p><textarea style="font-size: 17px" name="msg" rows="5" cols="50">
chr2:50843321-51318853
2p16.3(50843321-51318853)
2p16.3(50843321-51318853)x3 (advised format) </textarea></p></center><br>
		<br />
      </div>
      <script>
      function yesnoCheck(that){
       if (that.value == "breakpoint") {
         document.getElementById("ifmeh").style.display = "block";
         document.getElementById("ifins").style.display = "none";
       } else if (that.value == "cnv") {
         document.getElementById("ifmeh").style.display = "none";
         document.getElementById("ifins").style.display = "block";
       } else if (that.value == "bpcnv") {
         document.getElementById("ifmeh").style.display = "block";
         document.getElementById("ifins").style.display = "block";
       }
     }
</script> 
<p><center><input type="submit" value="Submit"/p></center><br>
<br />
</form>
<aaa><center>If you using this tool please acknowledge either by <i>This table was performed by the CNV-Content tool</i> or by citing <a href="http://www.insa.min-saude.pt/category/areas-de-atuacao/genetica-humana/">our reference publication</a></center><br><br />
<b><center><a href=/tadgctV2_tutorial.pdf>Reference manual</a></center></b><br><br />
<center><address>
Correspondance: <a href="mailto:doencasgenomicas@insa.min-saude.pt">Genomic Diseases Group</a>.<br><br />
</address>
<center><aaa><a href="http://www.insa.min-saude.pt/category/areas-de-atuacao/genetica-humana/">Department of Human Genetics</a></aaa></center><br>
<p>National Institute of Health Doutor Ricardo Jorge</p> </aaa></center>
<center><img src="https://cld.pt/dl/download/1f715328-21eb-49bd-b04c-b46bf2f08c61/aaa.jpg" width=641 height=122 border=0 alt=""><br><br />
<center><p><rodape>This file was last modified 13/12/2018</p></font></rodape>
</body>
</html>"""

showForm1()

