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

function highlightSG(){
	var i
	var j
	var seedR = []
	for (i=0; i< projects.responseJSON.projects[$("#txtPro").val()]["regions"].length; i++){
		for (j=0; j< projects.responseJSON.projects[$("#txtPro").val()]["regions"][i]["seed_region_names"].length; j++)
			if (seedR.includes(projects.responseJSON.projects[$("#txtPro").val()]["regions"][i]["seed_region_names"][j])==false){
				seedR.push(projects.responseJSON.projects[$("#txtPro").val()]["regions"][i]["seed_region_names"][j])
			}
	}
	var seedG = []
	for (i=0; i< seedR.length; i++){
		if (seedG.includes(projects.responseJSON.regions[seedR[i]]["seed_group"])==false){
			seedG.push(projects.responseJSON.regions[seedR[i]]["seed_group"])
		}


	}
	$("#txtSG").val(seedG)
}



function fillSG() {
	$('#txtSG')[0].innerHTML = ''
	var i;
	var grSelect = $('#txtSG')[0];
	var sCheck = []; //array to check if seed group is already added
	var regArr = []; //store regions in array
	var reg = projects.responseJSON.regions
	for (reg in projects.responseJSON.regions){
		regArr.push(reg)
	}
	for (i=0; i< regArr.length; i++){
		var option = document.createElement("option");	
		option.text = projects.responseJSON.regions[regArr[i]]["seed_group"];
		if (sCheck.includes(projects.responseJSON.regions[regArr[i]]["seed_group"]) == false && projects.responseJSON.regions[regArr[i]]["seed_group"] !== null){	
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



function deta(x){ // fills in the details of project/region
	var proj 
	proj = projects.responseJSON.projects[$('#txtPro').val()]
	
	var checkfillD
	if (checkfillD = 1){
		$('#txtName')[0].innerHTML ='';
		$('#txtDes')[0].innerHTML = '';
		$('#txtVar')[0].innerHTML = '';
	}
	if ($("#txtSR").val() == null || $("#txtCR").val() == null || $('#txtPro').val()== null){
		alert("No selection made.")
	}else if($("#txtPro").val().length>1){ 
		alert("More than one project selected.")

	}else {

		if (x == "CR"){
			if ($("#txtCR").val().length==1){
				var reg =  projects.responseJSON.regions[$('#txtCR').val() ]
				$("#txtSeq")[0].innerHTML = reg["reference"]
				$('#txtReg').val("coord")
			} else{
				alert("More than one coordinate selected or coordinate not selected.")
			}
		}else if (x == "SR"){
			if ($("#txtSR").val().length==1){
				var reg =  projects.responseJSON.regions[$('#txtSR').val()]
				$("#txtSeq")[0].innerHTML = reg["reference"]
				$('#txtReg').val("seed")
			} else{
				alert("More than one coordinate selected or coordinate not selected.")
			}
		}
		$('#txtName')[0].innerHTML = $('#txtPro').val() ;
		$('#txtDes')[0].innerHTML = proj["description"] 
		$('#txtVar')[0].innerHTML = proj["max_variants"] 
		$('#txtNuc').val(reg["is_nucleotide"].toString())
		
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
	if($("#txtPro").val().length>1){ 
		alert("More than one project selected.")
	}
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

function sgSel(){

	var i
	var j
	var k
	var m
	var sCheck = []; 
	var regArr = [];
	var regions = projects.responseJSON.regions

	if ($('#txtSR option').length == 0 && $('#txtCR option').length == 0){
		alert ("Project not selected.")
		
	}else{
		for (regions in projects.responseJSON.regions){
			regArr.push(regions)
		}
		regArr
		for (i=0; i< regArr.length; i++){ 
			if (projects.responseJSON.regions[regArr[i]]["seed_group"] == $("#txtSG").val()){ 
				sCheck.push ([regArr[i]]);
			} 
		}
		$("#txtSR").val(sCheck)
		regselCo()
		sCheck
		var por = projects.responseJSON.projects; 
		var proj = []	
		var seedR = []	
		for (por in projects.responseJSON.projects){ 
			proj.push(por)
		} 

		for (i=0; i< proj.length; i++){
			for (j=0; j< projects.responseJSON.projects[proj[i]]["regions"].length; j++){
				for (k=0; k< projects.responseJSON.projects[proj[i]]["regions"][j]["seed_region_names"].length; k++){
					for (m=0; m<sCheck.length; m++){
						if (projects.responseJSON.projects[proj[i]]["regions"][j]["seed_region_names"][k].includes(sCheck[m])){
							if(seedR.includes(proj[i])==false){
								seedR.push(proj[i])
							}
						}
					}
				}
			}
		}
		$("#txtPro").val("")
		$("#txtPro").val(seedR)
	}
}



/*
function openForm(x) {
    document.getElementById(x).style.display = "block";
}
function closeForm(x) {
    document.getElementById(x).style.display = "none";
}
*/

