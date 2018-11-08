function addPR(){
	var newPR = $('#orReg').val();
	var newV = $('#arr').val();
	if ($('#proreg').val() == "projects"){ 
		projects.responseJSON.projects[newPR] = newV;
		} 
	} else if ($('#proreg').val() == "regions"){ 
		
		projects.responseJSON.regions[newPR] = newV;
	}
}

function delPR(){
	delete projects.responseJSON.($('#proreg').val())[$('#txtF').val()]
}
