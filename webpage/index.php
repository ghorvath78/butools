<?php
header("Cache-Control: no-cache, must-revalidate"); // HTTP/1.1
header("Expires: Sat, 26 Jul 1997 05:00:00 GMT"); // Date in the past
?>

<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <!-- The above 3 meta tags *must* come first in the head; any other head content must come *after* these tags -->
    <meta name="description" content="">
    <meta name="author" content="">
    <link rel="icon" href="btlogo48.png">
    <title>BuTools V2.0</title>
    <!-- Bootstrap core CSS -->
    <link href="css/bootstrap.min.css" rel="stylesheet">
    <!-- Custom styles for this template -->
    <link href="butools.css" rel="stylesheet">
  </head>
  <body>

<?php
  for ($i = 1; $i <= 6; $i++)
    $liclass[$i] = "";
  $liclass[3] = "class=\"dropdown\"";
  if (!isset($_GET["page"]) || $_GET["page"]<1 || $_GET["page"]>40)
    $active = 1;
 else
    $active = $_GET['page'];
 if ($active>30)
    $liclass[3] = "class=\"dropdown active\"";
 else
    $liclass[$active] = "class=\"active\"";
?>



<div class="container">

    <nav class="navbar navbar-inverse">
      <div class="container-fluid">
        <div class="navbar-header">
          <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
            <span class="sr-only">Toggle navigation</span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a class="navbar-brand" href="index.php?page=1">BuTools</span></a>
        </div>
        <div id="navbar" class="collapse navbar-collapse">
          <ul class="nav navbar-nav">
            <li <?php echo $liclass[1];?>><a href="index.php?page=1">News</a></li>
            <li <?php echo $liclass[2];?>><a href="index.php?page=2">About</a></li>
            <li <?php echo $liclass[6];?>><a href="index.php?page=6">Download</a></li>
            <li <?php echo $liclass[3];?>>
                <a href="#" class="dropdown-toggle" data-toggle="dropdown" role="button" aria-haspopup="true" aria-expanded="false">Demo <span class="caret"></span></a>
                <ul class="dropdown-menu">
                  <li><a href="index.php?page=31">PH package</a></li>
                  <li><a href="index.php?page=32">MAP package</a></li>
                  <li><a href="index.php?page=33">Trace package</a></li>
                  <li><a href="index.php?page=34">Fitting package</a></li>
                  <li><a href="index.php?page=35">MAM package</a></li>
                  <li><a href="index.php?page=36">Queues package</a></li>
                </ul>
              </li>
            <li <?php echo $liclass[4];?>><a href="doc/index.html">Documentation</a></li>
            <li <?php echo $liclass[5];?>><a href="index.php?page=5">Contact</a></li>
          </ul>
        </div><!--/.nav-collapse -->
      </div>
    </nav>


<?php
  $pages[1] = "news.php";
  $pages[2] = "about.php";
  $pages[3] = "demo.php";
  $pages[4] = "doc.php";
  $pages[5] = "contact.php";
  $pages[6] = "download.php";

  if ($active>30)
    include "demo.php";
  else
    include $pages[$active];
?>


</div><!-- /.container -->



    <!-- Bootstrap core JavaScript
    ================================================== -->
    <!-- Placed at the end of the document so the pages load faster -->
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.3/jquery.min.js"></script>
    <script src="js/bootstrap.min.js"></script>
  </body>
</html>

