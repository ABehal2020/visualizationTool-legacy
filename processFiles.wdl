task processFileRun {
	command {

	}
	output {

	}
	runtime {
		docker: 'pyconnectwdl:sv2'
	}
}

workflow computeNum {
	call processFileRun
}