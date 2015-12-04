<?php

// this prevents images loading when they're being generated  
$HEADER = true;
while (true) 
{
 $filesizes = array(filesize("V_full.png"), filesize("V.png"), filesize("T_full.png"), filesize("T.png"), filesize("S.png"));

 if (!in_array(0, $filesizes)) 
 {
   $HEADER = false;
   break; 
 } 

if (HEADER) header('Location: index.php');

}

// produce a tailed log file of T.dat. Remember, file must be owned by www-data
system('echo "Datetime\t\tCase V.\t\tGold-plate V.\tFloor V.\tCase T.\t\tGold-plate T.\tFloor T.\tThermo. Active?\tThermo. Sensor?\tThermo. Temp?" > tmp_TTail.txt', $retval);
system('echo "---------------------------------------------------------------------------------------------------------------------------------------------------------------------" >> tmp_TTail.txt', $retval);
system('tail -60 T.dat | awk -F \' \' \'{print strftime("%d/%m/%y %H:%M:%S",$1)"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t\t"$9"\t\t"$10}\' >> tmp_TTail.txt', $retval); 

?>
<head>
<META HTTP-EQUIV="Refresh" CONTENT="15">
</head>
<body>
<table align="center">
 <!--<tr><td colspan="2" style="text-align: center; font-size: 25px;"><u>Voltages</u></td></tr>
 <tr><td><img src="V_full.png"/></td>
 <td><img src="V.png"/></td></tr>
 <tr><td>&nbsp;</td></tr>-->
 <!--<tr><td colspan="2" style="text-align: center; font-size: 25px;"><u>Temperatures</u></td></tr>
 <tr><td><img src="T_full.png"/></td>-->
 <td><img src="T.png"/></td></tr>
 <td><img src="S.png"/></td></tr>
</table>
</body>
<br/>
<!--<div style="text-align: center;"><a href="tmp_TTail.txt">Link to tailed T.dat file</a></div>-->
</html>
