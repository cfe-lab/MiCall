function download(filename, text) {

  var element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
  element.setAttribute('download', filename);
  
  document.body.appendChild(element);

  element.click();
}




function fillPro() {  //fills in the list of projects
	$('#txtPro')[0].innerHTML = ''
	var proregSelect = $("#txtPro")[0]; 
	var por = projects.responseJSON.projects; 
	for (por in projects.responseJSON.projects){ 
		var option = document.createElement("option"); 
		option.text = por; 
		proregSelect.add(option); 
	} 
}

function fillSG() {
	$('#txtSG')[0].innerHTML = ''
	var i;
	var grSelect = $('#txtSG')[0];
	var sCheck = []; //array to check if seed group is already added
	var regArr = []; //store regions in array
	var regions = projects.responseJSON.regions
	for (regions in projects.responseJSON.regions){
		regArr.push(regions)
	}
	for (i=0; i< regArr.length; i++){
		var option = document.createElement("option");	
		option.text = projects.responseJSON.regions[regArr[i]]["seed_group"];
		if (sCheck.includes(projects.responseJSON.regions[regArr[i]]["seed_group"]) == false){	
			sCheck.push (projects.responseJSON.regions[regArr[i]]["seed_group"]);
			grSelect.add(option);
		} 
	}
}


function fillSR(){ //fills seed regions
	$('#txtSR')[0].innerHTML = ''
	var i;
	var j;
	var check = [];
	var region = projects.responseJSON.projects[$('#txtPro').val()].regions
	var seregSelect = $("#txtSR")[0];
	for (i=0; i< region.length; i++){
		var pro = $('#txtPro').val();
		var seedR = region[i]["seed_region_names"]
		for (j=0; j<seedR.length; j++){	
			var option = document.createElement("option");	
			option.text = seedR[j];
			if (check.includes(seedR[j]) == false){	
				check.push (seedR[j]);
				seregSelect.add(option);
			}
		}
	}
}

function fillCR(){ //fills in list of coordinate regions of a project
	$('#txtCR')[0].innerHTML = ''
	var i;
	var region = projects.responseJSON.projects[$('#txtPro').val()].regions
	var coregSelect = $("#txtCR")[0];
	for (i=0; i<region.length; i++){
		var pro = $('#txtPro').val();
		var option = document.createElement("option");
		option.text = region[i]["coordinate_region"]
		coregSelect.add(option);
	}

}



function deta(){ // fills in the details of project/region
	var proj 
	proj = projects.responseJSON.projects[$('#txtPro').val()]
	
	var checkfillD
	if (checkfillD = 1){
		$('#nameP')[0].innerHTML ='';
		$('#desc')[0].innerHTML = '';
	}
	if ($("#txtSR").val() == null && $("#txtCR").val() == null && $('#txtPro').val()== null){
		alert("Region has not been selected.")
	} else {

		if ($("#txtCR").val().length == 1){
			var reg =  projects.responseJSON.regions[$('#txtCR').val() ]
			$("#seq")[0].innerHTML = reg["reference"]
		}else if ($("#txtSR").val().length == 1){
			var reg =  projects.responseJSON.regions[$('#txtSR').val()]
			$("#seq")[0].innerHTML = reg["reference"]
		}
		$('#nameP')[0].innerHTML = $('#txtPro').val() ;
		$('#desc')[0].innerHTML = proj["description"] + "\n" + "Max Variants: " + proj["max_variants"] + "\n"+ "Is Nucleotide: " + reg["is_nucleotide"]
	}
	checkfillD = 1
}

/*function addPR(){   // to be fixed, not relevant
	var newPR = $('#orReg').val();
	var newV = $('#arr').val();
	var newValue = {"thing":$('#orReg').val(), "value":$('#arr').val()}
	if ($('#proreg').val() == "projects"){
		projects.responseJSON.projects[newPR] = newValue;
		fill();
	} else if ($('#proreg').val() == "regions"){
		projects.responseJSON.regions[newPR] = newValue;
		fill();
	}
	alert(newPR + " has been added to " + $('#proreg').val() + ".");
}

function delPR(){
	if (confirm("Are you sure you want to delete?")) {
		if ($('#proreg').val() == "projects"){
			delete projects.responseJSON.projects[$('#txtPro').val()]
			fill();
		} else if ($('#proreg').val() == "regions"){
			delete projects.responseJSON.regions[$('#txtPro').val()]
			fill();
		}
	}
}*/

function coselReg(){
	var region = projects.responseJSON.projects[$('#txtPro').val()].regions;
	var coordArr = [];
	var i;
	for (i = 0; i < region.length; i++){
		coordArr.push (region[i]["coordinate_region"]);
	}
	var arrNum = coordArr.indexOf($("#txtCR").val()[0])
	$("#txtSR").val(region[arrNum]["seed_region_names"])
}

function regselCo(){
	var region = projects.responseJSON.projects[$('#txtPro').val()].regions;
	var check = []
	for (i = 0; i < region.length; i++){
		if (region[i]["seed_region_names"].includes($("#txtSR").val()[0])){
			check.push (region[i]["coordinate_region"])
		}
	}
	$("#txtCR").val(check)
}

/*function groupsel(){
	var i;
	var regArr = [];
	var reg = projects.responseJSON.regions; 
	var regArr= [] 
	for (reg in projects.responseJSON.regions){
		regArr.push(reg)
	} 
	for(i=0; regArr.length ; i++){
		projects.responseJSON.regions[regArr[i]]
	
	}
}


function openForm(x) {
    document.getElementById(x).style.display = "block";
}
function closeForm(x) {
    document.getElementById(x).style.display = "none";
}
*/

