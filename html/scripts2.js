//loads the json file
$(window).on('load',function(){
	$("#inputfile").val("")
	sprojects = $.getJSON('projects.json')
	document.getElementById('inputfile').addEventListener('change', handleFileSelect, false);
})

function load(){ //for when there has not been a json file chosen, uses default projects.json
	projects = sprojects.responseJSON
	clears();
	fillPro();
	fillSG();
	blankR();
	$("#loadbut").show()
	$("#basejson").hide()
}

function handleFileSelect(evt) {
    /*
    Handle file selected by browser dialog.
     */	
	reader = new FileReader()
	files = evt.target.files; // FileList object
	f = files[0];
	
	// Only process image files.
	if (!f.type.match('application/json')) {
        alert('Sorry, this script will only process JSON files.  ' +
            'You attempted to process a file of type: ' + f.type);
		return;
	}
	
	
	reader.readAsText(f);
	reader.onload = fileReadComplete;  // need some event handler here
}

function fileReadComplete(e) {
	projects = JSON.parse(reader.result)	
	$("#loadbut").show()
	$("#basejson").hide()
	clears();
	fillPro();
	fillSG();
	blankR();
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

function download() {
	var prepop = ""
	if ($("#inputfile")[0].files[0] == null){
			prepop = "download.json"
		}else{
			prepop = $("#inputfile")[0].files[0].name
	}
	if (projects !== undefined){
		
		
		var name = prompt ("Enter the file name for download." , prepop)
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

