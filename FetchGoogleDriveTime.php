<?php
// This script retrieves travel times from Google's Distance Matrix Service.
// The script takes an csv files as input, the file contains the list of X (latitude), Y (longitude) coordinates. 
// Travel times and distance are written to csv files as upper-triangle matrices. 
// Symmetric travel times are returned. To achieve this, the smallest of the two travel times A->B and B->A is used.

	function curl_request($sURL,$sQueryString=null) {
		$cURL=curl_init();
		curl_setopt($cURL,CURLOPT_SSL_VERIFYPEER, false);
		curl_setopt($cURL,CURLOPT_URL,$sURL.'?'.$sQueryString);
		curl_setopt($cURL,CURLOPT_RETURNTRANSFER, TRUE);
		$cResponse=trim(curl_exec($cURL));		
		curl_close($cURL);
		return $cResponse;
	}
	
	function getDistance($x1, $y1, $x2, $y2) {
		$sResponse=curl_request('https://maps.googleapis.com/maps/api/distancematrix/json', 'origins='.$x1.','.$y1.'&destinations='.$x2.','.$y2.'&mode=driving&sensor=false&key=AIzaSyA2xXJSTRZY4h2zWD8GdfSrCUBGFq9AQ7o');		
		$oJSON=json_decode($sResponse);			
		//echo $oJSON->status;
		if ($oJSON->status=='OK') {
			$fDistance=(float)preg_replace('/[^\d\.]/','',$oJSON->rows[0]->elements[0]->distance->value);
			$fTime=(float)preg_replace('/[^\d\.]/','',$oJSON->rows[0]->elements[0]->duration->value);
		} else {
			$fDistance=0;
			$fTime=0;
		}
		//echo 'Distance: '.$fDistance.'<br />Duration: '.$fTime.PHP_EOL;
		return array($fDistance, $fTime);
	}
	
	function csvToArray($file='', $delimiter=',') {
		if (($handle = fopen($file, "r")) !== FALSE) { 
            $i = 0; 
            while (($lineArray = fgetcsv($handle, 4000, $delimiter)) !== FALSE) { 
                for ($j=0; $j<count($lineArray); $j++) { 
                    $data2DArray[$i][$j] = $lineArray[$j]; 
                } 
                $i++; 
            } 
            fclose($handle); 
        } 
        return $data2DArray;
	}
	
	function arrayToCsv($file='file.csv', $list) {
		$fp = fopen($file, 'w');
		foreach ($list as $fields) {
			fputcsv($fp, $fields);
		}
		fclose($fp);
	}	
	
	// HERE THE ACTUAL CODE BEGINS
	ini_set('max_execution_time', 1800);
	$sleepingTime = 3600*24+1;
	$limit = 2000;
	$count = 0;
	$coordinates = csvToArray('Locations.csv', ',');
		
	$length = sizeof($coordinates);
	
	for ($i=0; $i < $length; $i++) {
		for ($j=0; $j < $length; $j++) {
			$distance[$i][$j] = 0;
			$duration[$i][$j] = 0;
		}
	}
	
	for ($i=0; $i < $length; $i++) {	
		for ($j=$i+1; $j < $length; $j++) {				
			$tmp = getDistance($coordinates[$i][0], $coordinates[$i][1], $coordinates[$j][0], $coordinates[$j][1]);
			$tmp2 = getDistance($coordinates[$j][0], $coordinates[$j][1], $coordinates[$i][0], $coordinates[$i][1]);
			
			$distance[$i][$j] = min($tmp[0], $tmp2[0]);
			$duration[$i][$j] = min($tmp[1], $tmp2[1]);
			echo("[".$i.", ".$j."]");
			$count++;
			if ((($count % $limit) == 0) && ($count > 0)) {
				arrayToCsv('distance.csv', $distance);
				arrayToCsv('duration.csv', $duration);
				sleep($sleepingTime);
			}
		}
	}
	
	arrayToCsv('distance.csv', $distance);
	arrayToCsv('duration.csv', $duration);
	
?>