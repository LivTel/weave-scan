<?php

// this prevents images loading when they're being generated  
$HEADER = true;
while (true) 
{
 $filesizes = array(filesize("T.png"), filesize("D.png"), filesize("H.png"));

 if (!in_array(0, $filesizes)) 
 {
   $HEADER = false;
   break; 
 } 

if (HEADER) header('Location: index.php');

}

?>
<head>
<META HTTP-EQUIV="Refresh" CONTENT="15">
</head>
<body>
<table align="center">
 <tr><td><img style="width:480px;" src="T.png"/></td>
 <td><img style="width:480px;" src="H.png"/></td></tr>
 <td><img style="width:480px;" src="D.png"/></td></tr>
</table>
</body>
<br/>
</html>
