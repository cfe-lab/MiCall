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


function download() {
	if (projects !== undefined){
		name = prompt("Enter file name for download.")
		if(name !== "null"){
			if(name == ""){
				if ($("#inputfile")[0].files[0] == null){
					name = "download"
				}else{
					name = $("#inputfile")[0].files[0].name
				}
			}		
			if(name.indexOf(".json") == -1){
				name += ".json"
			}
			var element = document.createElement('a');
			element.setAttribute('href', 'data:text/plain;charset=utf-8,' + encodeURIComponent(JSON.stringify(projects,null,"\t")));
			element.setAttribute('download', name);
		  	  element.setAttribute("target", "_blank")
			document.body.appendChild(element);

			element.click();
		}
	}else{
		alert("No data to be downloaded.")
	}
}
