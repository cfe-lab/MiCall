function launchSpec(dataProvider)
{
    var ret = {
        commandLine: [ "basespace" ],
        // NOTE: Replace the tag below with the specific release version (e.g. :v4.10.0)
        // before creating a new BaseSpace app version.  Illumina recommends pinning
        // Native Apps to a specific image digest rather than a mutable tag for
        // reproducibility (see Illumina Native App publishing guide).
        // Example digest format: "docker.illumina.com/cfe_lab/micall@sha256:<hash>"
        containerImageId: "docker.illumina.com/cfe_lab/micall:vX.Y.Z",
        Options: [ "bsfs.enabled=true" ]
    };
    return ret;
}

// example multi-node launch spec
/*
function launchSpec(dataProvider)
{
    var ret = {
        nodes: []
    };
    
    ret.nodes.push({
        appSessionName: "Hello World 1",
        commandLine: [ "cat", "/illumina.txt" ],
        containerImageId: "basespace/demo",
        Options: [ "bsfs.enabled=true" ]
    });
    
    ret.nodes.push({
        appSessionName: "Hello World 2",
        commandLine: [ "cat", "/illumina.txt" ],
        containerImageId: "basespace/demo",
        Options: [ "bsfs.enabled=true" ]
    });
    
    return ret;
}
*/

/* 
function billingSpec(dataProvider) {
    return [
    {
        "Id" : "insert product ID here",
        "Quantity": 1.0
    }];
}
*/
