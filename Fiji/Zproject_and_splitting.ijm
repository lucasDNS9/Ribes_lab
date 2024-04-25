//Elie BALLOUL 2015 elie.balloul@logoss.eu

//In the middle part of the file you can change the commands
//And put new ones from the record macro command

IGet=getDirectory("Which directory to treat?");
sep=File.separator();
ISave=File.getParent(IGet) + sep + "split_channels" + File.getName(IGet) + sep; 
if(File.isDirectory(ISave)==0) File.makeDirectory(ISave);

function FindAndSave(Input, Output) {
	ListFile=getFileList(Input);
	for(n=0; n<ListFile.length; n++) {
		Target=Input + ListFile[n];
		if(File.isDirectory(Target)) {
			ITSave=Output + sep + ListFile[n] + sep; 
			if(File.isDirectory(ITSave)==0) File.makeDirectory(ITSave);
			FindAndSave(Target, ITSave);	
		} else {
			open(Target);
			varName=getTitle();
			OriCurrent=getImageID();

			//////////////
			///Change the next commands according to your markers
			//////////////
			title = getTitle();
			run("Z Project...", "projection=[Average Intensity]");
			run("Split Channels");

			//Channel 1
			selectWindow("C1-AVG_"+title);
			rename("C1");
			run("Red");
			saveAs("Tiff", ISave + "C1-PH3-"+title);
			//waitForUser("set the color then click on OK");

			//Channel 2 
			selectWindow("C2-AVG_"+title);
			rename("C2");
			run("Green");
			saveAs("Tiff", ISave + "C2-PAX3-"+title);
			//waitForUser("set the color then click on OK");

			//Channel 3
			selectWindow("C3-AVG_"+title);
			rename("C3");
			run("Blue");
			saveAs("Tiff", ISave + "C3-DAPI-"+title);
			//waitForUser("set the color then click on OK");

			//Channel 4
			//selectWindow("C4-AVG_"+title);
			//rename("C4");
			//run("Blue");
			//saveAs("Tiff", ISave + "C4-Dapi-"+title);
			//waitForUser("set the color then click on OK");

			/////// The end

			close("*");
		}
	}
}

FindAndSave(IGet, ISave);
