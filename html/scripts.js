//loads the json file
projects = $.getJSON('projects.json')


function loadChange(){
	$("#loadbut")[0].innerHTML = "Refresh"
}

function download(filename, text) { //downloads the json file
  var element = document.createElement('a');
  element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(text));
  element.setAttribute('download', filename);
  
  document.body.appendChild(element);

  element.click();
}


function fillPro() {  //fills in the list of projects
	$('#txtPro')[0].innerHTML = ''
	var proSelect = $("#txtPro")[0]; 
	var por = projects.responseJSON.projects; 
	for (por in projects.responseJSON.projects){ 
		var option = document.createElement("option"); 
		option.text = por; 
		proSelect.add(option); 
	} 
}


function highlightSG(){ //highlights related seed groups when project is chosen
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



function fillSG() { //fills list of seed groups
	$('#txtSG')[0].innerHTML = ''
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
		if (sCheck.includes(projects.responseJSON.regions[regArr[i]]["seed_group"]) == false && projects.responseJSON.regions[regArr[i]]["seed_group"] !== null && projects.responseJSON.regions[regArr[i]]["seed_group"] !== "null"){	
			sCheck.push (projects.responseJSON.regions[regArr[i]]["seed_group"]);
			grSelect.add(option);
		} 
	}
}

function blankR(){ //clears both region select fields
	$("#txtCR")[0].innerHTML =""
	$("#txtSR")[0].innerHTML =""
}

function fillSR(){ //fills seed regions
	$('#txtSR')[0].innerHTML = ''
	var check = [];
	var seregSelect = $("#txtSR")[0];
	if ($("#txtPro").val().length ==1){
		var region = projects.responseJSON.projects[$('#txtPro').val()].regions
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
	}else if ($("#txtPro").val().length >1){
		for (k=0; k<$("#txtPro").val().length ; k++){
		var region = projects.responseJSON.projects[$('#txtPro').val()[k]].regions
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
	}
}

function fillCR(){ //fills in list of coordinate regions of a project
	$('#txtCR')[0].innerHTML = ''
	var coregSelect = $("#txtCR")[0];
	var check = [];
	if ($("#txtPro").val().length ==1){
		var region = projects.responseJSON.projects[$('#txtPro').val()].regions
		for (i=0; i<region.length; i++){
			var pro = $('#txtPro').val();
			var option = document.createElement("option");
			option.text = region[i]["coordinate_region"]
			coregSelect.add(option);
		}
	}else if ($("#txtPro").val().length >1){
		for (k=0; k<$("#txtPro").val().length ; k++){
		var region = projects.responseJSON.projects[$('#txtPro').val()[k]].regions
			for (i=0; i<region.length; i++){
				var pro = $('#txtPro').val();
				var option = document.createElement("option");
				option.text = region[i]["coordinate_region"]
				if (check.includes(region[i]["coordinate_region"]) == false){	
						check.push (region[i]["coordinate_region"]);
						coregSelect.add(option);
				}
			}
		}
	}
}


function coselReg(){ //highlights related seed region fields when coordinate region is selected
	if ($("#txtPro").val().length ==1){
		var region = projects.responseJSON.projects[$('#txtPro').val()].regions;
		var coordArr = [];
		for (i = 0; i < region.length; i++){
			coordArr.push (region[i]["coordinate_region"]);
		}
		var arrNum = coordArr.indexOf($("#txtCR").val()[0])
		var seedRG = region[arrNum]["seed_region_names"]
		$("#txtSR").val(seedRG)
		if (region[arrNum]["seed_region_names"].length !== 0){
			$("#txtSG").val(projects.responseJSON.regions[seedRG[0]]["seed_group"])
		}else{
			$("#txtSG").val("")
		}
	}else if ($("#txtPro").val().length >1){
		var check = []
		for (j=0; j<$("#txtPro").val().length ; j++){
			var region = projects.responseJSON.projects[$('#txtPro').val()[j]].regions;
			var coordArr = [];
			for (i = 0; i < region.length; i++){
				coordArr.push (region[i]["coordinate_region"]);
			}
			var arrNum = coordArr.indexOf($("#txtCR").val()[0])
			if (arrNum !== -1){
				for (k=0; k<region[arrNum]["seed_region_names"].length; k++){
					if (check.includes(region[arrNum]["seed_region_names"][k].toString()) == false){
						check.push(region[arrNum]["seed_region_names"][k].toString())
					}	
				}		
			}	
		}
		$("#txtSR").val(check)
		$("#txtSG").val(projects.responseJSON.regions[check[0]]["seed_group"])
	}
}

function regselCo(){ //highlights related coordinate region fields when seed region is selected
	if ($("#txtPro").val().length ==1){
		var region = projects.responseJSON.projects[$('#txtPro').val()].regions;
		var check = []
		for (i = 0; i < region.length; i++){
			if (region[i]["seed_region_names"].includes($("#txtSR").val()[0])){
				check.push (region[i]["coordinate_region"])
			}
		
		}
	$("#txtSG").val(projects.responseJSON.regions[$("#txtSR").val()]["seed_group"])		
	}else if ($("#txtPro").val().length >1){
		var check = []
		for (j=0; j<$("#txtPro").val().length ; j++){
			var region = projects.responseJSON.projects[$("#txtPro").val()[j]].regions;
			
			for (i = 0; i < region.length; i++){
				if (region[i]["seed_region_names"].includes($("#txtSR").val()[0])){
					check.push (region[i]["coordinate_region"])
				}
			}
		}
	 $("#txtSG").val(projects.responseJSON.regions[$("#txtSR").val()[0]]["seed_group"])
	}

	$("#txtCR").val(check)
}

function sgSel(){ //select seed group, finds and highlights projects and regions related
	var sCheck = []; 
	var regArr = [];
	var regions = projects.responseJSON.regions

		for (regions in projects.responseJSON.regions){
			regArr.push(regions)
		}
		for (i=0; i< regArr.length; i++){ 
			if (projects.responseJSON.regions[regArr[i]]["seed_group"] == $("#txtSG").val()){ 
				sCheck.push ([regArr[i]]);
			} 
		}
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
		fillSR()
		fillCR()
		$("#txtSR").val(sCheck)
		regselCo()
		
}

function deta(x){ // fills in the details of project/region
	var proj 
	proj = projects.responseJSON.projects[$('#txtPro').val()]
	$('#txtNuc').val("");
	
	if ($("#txtSR").val() == null || $("#txtCR").val() == null || $('#txtPro').val()== null){
		alert("No selection made.")
	}else {

		if (x == "CR"){
			if ($("#txtCR").val().length==1){
				var reg =  projects.responseJSON.regions[$('#txtCR').val() ]
				$("#txtSeq").val(reg["reference"].toString().replace(/\,/g, ''))
				$('#txtReg').val("coordinate")
				$('#txtRN').val($("#txtCR").val())
				regType = "coordinate" //global variable to check region type
			} else{
				alert("More than one coordinate selected or coordinate not selected.")
			}
		}else if (x == "SR"){
			if ($("#txtSR").val().length==1){
				var reg =  projects.responseJSON.regions[$('#txtSR').val()]
				$("#txtSeq").val(reg["reference"].toString().replace(/\,/g, ''))
				$('#txtReg').val("seed")
				$('#txtRN').val($("#txtSR").val())
				regType = "seed" //global variable to check region type
			} else{
				alert("More than one coordinate selected or coordinate not selected.")
			}
		}
		globReg = reg //global variable to store region chosen
		$('#txtNuc').val(reg["is_nucleotide"].toString())
		if (reg["is_nucleotide"] == false){
			$("#txtSeG").val("null")
		}else{
			$("#txtSeG").val(reg["seed_group"])
		}
		
	}
	if($("#txtPro").val().length == 1){
		$("#txtName").val($('#txtPro').val()) ;
		$('#txtDes').val(proj["description"]) 
		$('#txtVar').val(proj["max_variants"]) 

	} else{
		$('#txtName').val( "Project not selected.")
		$('#txtDes').val("N/A")
		$('#txtVar').val("N/A")

	}
}

function clears(){ //clears the details fields
	$("#txtName").val("")
	$('#txtDes').val("")
	$('#txtVar').val("")
	$("#txtSeG").val("")
	$('#txtNuc').val("")
	$("#txtSeq").val("")
	$('#txtReg').val("")
	$('#txtRN').val("")
}


function editCheck(){ //check to see which fields are edited and if changes are permitted
	var proj = projects.responseJSON.projects[$('#txtPro').val()]
	var reg =  projects.responseJSON.regions
	var changes =""
	var matchOK 
	

	if ($("#txtPro").val().length > 1 || $("#txtPro").val().length == 0){ //check for name change
		alert ("Project needs to be selected.")
		matchOK = false
	}else if($("#txtSR").val() == "" || $("#txtSR").val() == ""){
		alert("Region needs to be selected.")
	}else if($("#txtPro").val().toString() == $("#txtName").val()){
		newProj = $("#txtPro").val().toString()
	}else if($("#txtPro").val().toString() !== $("#txtName").val()){
		newProj = $("#txtName").val()
		changes += "\n" + "Project Name"
	}

	
	if ($("#txtVar").val() == proj["max_variants"]){ //check for max variant change
		newVar = proj["max_variants"]
	}else if ($("#txtVar").val() !== proj["max_variants"]){	
		newVar = $("#txtVar").val()
		changes += "\n" + "Max Variants"
	}
	

	if ($("#txtDes").val() == proj["description"]){ //check for description change
		newDes = proj["description"]
	}else if ($("#txtDes").val() !== proj["description"]){
		newDes = $("#txtDes").val()
		changes += "\n" + "Description"
	}

	
	if ($("#txtReg").val() == regType){ //check for region type
		newregType = regType
	}else if ($("#txtReg").val() !== regType){
		newregType = $("#txtReg").val()
		changes += "\n" + "Region Type"
	}


	if (regType == "coordinate"){ //check for region name
		regName = $("#txtCR").val().toString()
	}else if (regType == "seed"){
		regName = $("#txtSR").val().toString()
	}
	if ($("#txtRN").val() !== regName){
		changes += "\n" + "Region Name"
	} 
 	newregName = $("#txtRN").val()


	//check for seed group
	if ($("#txtReg").val() == "coordinate" && $("#txtSeG").val() !== "null"){
		alert("seed group for coordinate regions is 'null'")
		$("#txtSeG").val("null")
		newSG = $("#txtSeG").val()
	}else if ($("#txtReg").val() == "coordinate" && $("#txtSeG").val() == "null"){
		newSG = "null"
	}else if ($("#txtReg").val() == "seed" && $("#txtSeG").val() == globReg["seed_group"]){
		newSG = $("#txtSeG").val()
	}else if ($("#txtReg").val() == "seed" && $("#txtSeG").val() !== globReg["seed_group"]){
		newSG = $("#txtSeG").val()
		changes += "\n" + "Seed Group"
	}



	if (globReg["reference"].toString().replace(/\,/g, '') == $("#txtSeq").val().replace(/\,/g, '').toUpperCase()){ //check for sequence
		newSeq = globReg["reference"]
	}else if (globReg["reference"].toString().replace(/\,/g, '') !== $("#txtSeq").val().replace(/\,/g, '').toUpperCase()){ //converts sequence to array
		newSeq = []
		var arrLength = Math.ceil($("#txtSeq").val().replace(/\,/g, '').length/65) 
		for (i=0; i< arrLength-1; i++){
			var seqBlock = ""
			for (j=0; j<65; j++){
				if ($("#txtSeq").val()[i*65 +j].toUpperCase().length!== undefined){
					seqBlock += $("#txtSeq").val()[i*65 +j].toUpperCase()
				}
			}
			newSeq[i] = seqBlock
		}	
		var seqBlock = ""
		for (i=(arrLength-1)*65; i<$("#txtSeq").val().replace(/\,/g, '').length; i++){
			seqBlock += $("#txtSeq").val()[i].toUpperCase()
		}
		newSeq[arrLength-1] = seqBlock
			
		changes += "\n" + "Reference Sequence"
	}
	isNuc = 0
	
	for (i=0; i<newSeq.toString().length; i++){ //checks if sequence is nucleotide or not
		if(newSeq.toString()[i] !== "A" && newSeq.toString()[i] !== "T" && newSeq.toString()[i] !== "G" && newSeq.toString()[i] !== "C" && newSeq.toString()[i] !== "," && newSeq.toString()[i] !== "*" ){
			isNuc += 1
		}
	}
	
	//check if sequence matches region type
	if ($("#txtReg").val() == "coordinate"){
		if($("#txtNuc").val() == "false" && isNuc >0){
			matchOK = true

		}else{
			alert("For a coordinate region, sequence must be an Amino Acid Sequence")
			matchOK = false
		}
	}else if ($("#txtReg").val() == "seed"){
		if($("#txtNuc").val() == "true" && isNuc == 0){
			matchOK = true
		}else{
			alert("For a seed region, sequence must be a Nucleotide Sequence")
			matchOK = false
		}
	}
	if (isNaN($("#txtVar").val())){
		alert("Max Variants must be a number")
		matchOK = false
	}else{		
 	}
	if (matchOK == true && ($("#txtSR").val() !== "" || $("#txtSR").val() !== "")){
		if (changes == ""){
			if (confirm ("There were no changes")){
			}	
		}else{ 
			if (confirm ("The following will be changed:" + changes)){
				editChange()
			}	
		}
	}
}

function editChange(){ //makes the edits
	var proje = projects.responseJSON.projects
	var regi = projects.responseJSON.regions
	
	if (newProj !== $("#txtPro").val().toString()){
		proje[newProj] = proje[$("#txtPro").val().toString()]
		delete proje[$("#txtPro").val().toString()]
	}
	proje[newProj]["max_variants"] = newVar
	proje[newProj]["description"] = newDes
	

	
	var regionIndex = []
	var seedIndex
	var coorIndex = []
	
	if (newregType == "coordinate"){ //for when coordinate name is changed
		regi[newregName] = {"is_nucleotide": false, "reference": newSeq, "seed_group": newSG}
		for (i=0; i<$("#txtCR option").length -1;i++){		
			if ($("#txtCR option")[i].innerHTML == $("#txtCR").val()){
				coorIndex.push(i)
			}
		}
		proje[newProj]["regions"][coorIndex[0]]["coordinate_region"] = newregName

	}else if (newregType == "seed"){ //for when seed name is changed 
		regi[newregName] = {"is_nucleotide": true, "reference": newSeq, "seed_group": newSG}
		for (i=0; i<$("#txtCR").val().length; i++){
			for (j=0; j<$("#txtCR option").length -1;j++){		
				if ($("#txtCR option")[j].innerHTML == $("#txtCR").val()[i]){
					regionIndex.push(j)
				}
			}
			seedIndex = proje[newProj]["regions"][regionIndex[i]]["seed_region_names"].indexOf(regName)
			proje[newProj]["regions"][regionIndex[i]]["seed_region_names"][seedIndex] = newregName
		}	
	}

}

function delPR(){  //to delete a project or region
	var por = projects.responseJSON.projects;
	var check = []
	if ($("#txtPro").val().length !== 1){
		alert ("Project must be selected.")
	}else{
		if (confirm("Are you sure you want to delete?")) {
			if ($('#delType').val() == "project"){
				delete projects.responseJSON.projects[$('#txtPro').val()]
				alert ("project has been deleted")
			}else if ($('#delType').val() == "region"){ 
				var coorIndex = []
				var upPro = projects.responseJSON.projects[$('#txtPro').val()]["regions"]
				if (regType == "coordinate"){ 
					for (i=0; i<$("#txtCR option").length;i++){		
						if ($("#txtCR option")[i].innerHTML == $("#txtCR").val()){
							coorIndex.push(i)
						}
					}
					delete projects.responseJSON.projects[$('#txtPro').val()]["regions"][coorIndex[0]]
					
					//shift array indexes down
					for (i = coorIndex[0]; i<upPro.length; i++){
						upPro[i] = upPro[i+1]
					}
					upPro.length = upPro.length - 1
					for (por in projects.responseJSON.projects){
						for(i=0; i<projects.responseJSON.projects[por]["regions"].length; i++){
							if (projects.responseJSON.projects[por]["regions"][i]["coordinate_region"]==$("#txtCR").val()){
								check.push($("#txtCR").val())
							}
						}
					}
					if (check.length == 0){
						delete projects.responseJSON.regions[$("#txtCR").val()]
					}
				}else if (regType == "seed"){ 
					var regionIndex = [] //where in the project regions it is
					var seedIndex //index of where the seed region is in the specific region
					for (i=0; i<$("#txtCR").val().length; i++){
						for (j=0; j<$("#txtCR option").length;j++){		
							if ($("#txtCR option")[j].innerHTML == $("#txtCR").val()[i]){
								regionIndex.push(j)
							}
						}
						seedIndex = upPro[regionIndex[i]]["seed_region_names"].indexOf($("#txtSR").val().toString())
						delete upPro[regionIndex[i]]["seed_region_names"][seedIndex]
						for (k = seedIndex; k < upPro[regionIndex[i]]["seed_region_names"].length; k++){
							upPro[regionIndex[i]]["seed_region_names"][k] = upPro[regionIndex[i]]["seed_region_names"][k+1]
						}
					upPro[regionIndex[i]]["seed_region_names"].length = upPro[regionIndex[i]]["seed_region_names"].length -1
					}
					for (por in projects.responseJSON.projects){
						for(i=0; i<projects.responseJSON.projects[por]["regions"].length; i++){
							for(j=0; j<projects.responseJSON.projects[por]["regions"][i]["seed_region_names"].length; j++){
								if (projects.responseJSON.projects[por]["regions"][i]["seed_region_names"][j]==$("#txtSR").val()){
									check.push($("#txtSR").val())
								}
							}
						}
					}
					if (check.length == 0){
						delete projects.responseJSON.regions[$("#txtSR").val()]
					}	
				}
				
				

				alert (regType + " region has been deleted.")
			}
			clears()
			fillPro()
			fillSG()
			blankR()
		}
	}
}

function fillAddPro() {  //fills in the list of projects
	$('#addProj')[0].innerHTML = ""
	var proSelect = $('#addProj')[0]; 
	var por = projects.responseJSON.projects; 
	var option = document.createElement("option"); 
		option.text = "New Project";
		proSelect.add(option)
	for (por in projects.responseJSON.projects){ 
		var option = document.createElement("option"); 
		option.text = por; 
		proSelect.add(option); 
	} 
}

function fillAddCR(){ //fills coordinate regions in add form
	$("#cooregion")[0].innerHTML = ""
	var pro = projects.responseJSON.projects
	if ($("#addProj").val() !== "New Project"){
		var cooSelect = $("#cooregion")[0]
		for (i=0; i< pro[$("#addProj").val()]["regions"].length; i++){
			var option = document.createElement("option"); 
			option.text = pro[$("#addProj").val()]["regions"][i]["coordinate_region"]; 
			cooSelect.add(option); 
		} 
	}
}

function addPCheck(){ //checks for new or existing project
	
	var pro = projects.responseJSON.projects
	var reg = projects.responseJSON.regions
	if ($("#addProj").val() !== "New Project"){ //existing, fields cant be edited, fills the fields with existing data
		$("#addProjName").val($("#addProj").val())
		$("#addProjName").prop("readonly", true)
		$("#addVar").val(pro[$("#addProj").val()]["max_variants"].toString())
		$("#addVar").prop("readonly", true)
		$("#addDes").val(pro[$("#addProj").val()]["description"].toString())
		$("#addDes").prop("readonly", true)	

	}else{ //new project, fields are editable
		$("#addProjName").val("")
		$("#addProjName").prop("readonly", false)
		$("#addVar").val("")
		$("#addVar").prop("readonly", false)
		$("#addDes").val("")
		$("#addDes").prop("readonly", false)
	}
}

function addPO(){ //regarding the project fields
	if ($("#addProj").val() == "New Project"){
		var newProj = $("#addProjName").val()
		var newValue = {"max_variants": $("#addVar").val(), "description": $("#addDes").val(), "regions": []}
		var pro = projects.responseJSON.projects

		//checks for empty fields, and if max variants is a number
		if ($("#addProjName").val() !== "" && isNaN($("#addVar").val())== false && $("#addVar").val() !== "" && $("#addDes").val() !==""){ 
			pro[newProj] = newValue
			$("#addregion").show()
			$("#addProjName").prop("readonly", true)
			$("#addVar").prop("readonly", true)
			$("#addDes").prop("readonly", true)
		}else if(isNaN($("#addVar").val())== true){
			alert("Max Variants must be a number")
		}else if($("#addProjName").val() == "" || $("#addVar").val() == "" || $("#addDes").val() ==""){
			alert("Fields are empty.")
		}
		
	}else if ($("#addProj").val() !== "New Project"){
		$("#addregion").show()
	}
	
}


function addRCheck(){ //checks for coordinate or seed that user wanted to add
	if($("#addRegi").val() == "coordinate"){
		$("#seedgrouptext").hide()
	}else if ($("#addRegi").val() == "seed"){
		$("#seedgrouptext").show()
	}
	$("#common").show()
	$("#regionbut").show()
	$("#addRegName").val("")
	$("#addSequence").val("")
	$("#addSeedGroup").val("")
	
}


function seedShow(){ //displays the seed regions of the selected coordinate region in add form
	$("#seeregion")[0].innerHTML = ""
	var pro = projects.responseJSON.projects
	if ($("#addProj").val() !== "New Project"){
		var seeSelect = $("#seeregion")[0]
		for (i=0; i< pro[$("#addProj").val()]["regions"].length; i++){
			if (pro[$("#addProj").val()]["regions"][i]["coordinate_region"] == $("#cooregion").val()){		
				for (j=0; j< pro[$("#addProj").val()]["regions"][i]["seed_region_names"].length; j++){
					var option = document.createElement("option"); 
					option.text = pro[$("#addProj").val()]["regions"][i]["seed_region_names"][j]; 
					seeSelect.add(option); 
				}
			}
		} 
	}
}

function addRe(){
	var pro = projects.responseJSON.projects
	var reg = projects.responseJSON.regions
	var isNuc = 0	
	var exist = 0
	var addSeq = $("#addSequence").val().replace(/\,/g, '').toUpperCase()
	for (i=0; i<addSeq.length; i++){ //checks if sequence is nucleotide or not
		if(addSeq[i] !== "A" && addSeq[i] !== "T" && addSeq[i] !== "G" && addSeq[i] !== "C" && addSeq[i] !== "," && addSeq[i] !== "*" ){
			isNuc += 1  
		}
	}
	if($("#addRegi").val() == "coordinate"){
		for(i=0; i< $("#cooregion option").length; i++){
			if ($("#addRegName").val() == $("#cooregion option")[i].innerHTML){
				exist += 1
			}
		}
		if ($("#addRegName").val() !== "" && $("#addSequence").val() !== "" && exist == 0){//add check for existing coord with same name in projects and regions
			if(isNuc !== 0){
				pro[$("#addProjName").val()]["regions"].push ({"coordinate_region": $("#addRegName").val() ,"seed_region_names":""})
				reg[$("#addProjName").val()] = {"is_nucleotide": false, "reference": $("#addSequence").val(), "seed_group": null}
				fillAddCR()			
			}else{
				alert("Must be an amino acid sequence for coordinate region.")
			}
				
		}else if ($("#addRegName").val() == "" && $("#addSequence").val() == ""){
			alert("Fields must be filled in.")
		}else if (exist > 0){
			alert("Region already exists.")
		}
	}else if ($("#addRegi").val() == "seed"){ //add check for existing seed with same name in projects and regions
		var proIndex = []
		var seedStore = []
		for(i=0; i<$("#cooregion option:selected").length; i++){ //gets index of the chosen coordinate regions in project
			proIndex.push($("#cooregion option:selected")[i].index)
			for(j=0; j<pro[$("#addProjName").val()]["regions"][proIndex[i]]["seed_region_names"].length; j++){
				seedStore.push(pro[$("#addProjName").val()]["regions"][proIndex[i]]["seed_region_names"][j])
			}
		}
		for(i=0; i<seedStore.length; i++){
			if ($("#addRegName").val() == seedStore[i]){
				exist += 1
			}
		}

		var cooCheck = []
		if ($("#addRegName").val() !== "" && $("#addSequence").val() !== "" && $("#addSeedGroup").val() !== "" && exist == 0){
			for (i=0; i< $("#cooregion").val().length; i++){
				for (j=0; j<pro[$("#addProjName").val()]["regions"].length ; j++){ 
					if ($("#cooregion").val()[i] == pro[$("#addProjName").val()]["regions"][j]["coordinate_region"]){
						cooCheck.push(j)
					}	
				}
			}
			if(isNuc == 0){
				for (i=0; i< cooCheck.length; i++){
					pro[$("#addProjName").val()]["regions"][cooCheck[i]]["seed_region_names"].push($("#addRegName").val());
					reg[$("#addRegName").val()] = {"is_nucleotide": true, "reference": $("#addSequence").val(), "seed_group": $("#addSeedGroup").val()};
					console.log(projects.responseJSON.regions[$("#addProjName").val()])
				}
			}else{
				alert("Must be a nucleic acid sequence for seed region.")
			}
		}else if ($("#addRegName").val() == "" && $("#addSequence").val() == "" && $("#addSeedGroup").val() == ""){
			alert("Fields must be filled in.")
		}else if (exist > 0){
			alert("Region already exists.")
		}
	}

}

function openForm(x) { //open add or delete form
	document.getElementById(x).style.display = "block";
	$("#addProjName").val("")
	$("#addVar").val("")
	$("#addDes").val("")
	$("#addRegi").val("")
}
function closeForm(x) { //close add or delete form
	document.getElementById(x).style.display = "none";
}


