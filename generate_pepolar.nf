nextflow.preview.dsl = 2

usage = file("${workflow.scriptFile.getParent()}/usage/generate_pepolar")
engine = new groovy.text.SimpleTemplateEngine()

bindings = [
            "bold_tags": params.bold_tags,
            "sbref_tag": params.sbref_tag,
            "subjects": params.subjects,
            "rewrite": params.rewrite,
            "usebold": params.usebold
           ]

toprint = engine.createTemplate(usage.text).make(bindings)
printhelp = params.help

// Basic argument check
req_param = ["--study": params.study,
             "--bold_tags": params.bold_tags,
             "--out": params.out]
missing_arg = req_param.grep{ (it.value == null || it.value == "") }
if (missing_arg){
    log.error("Missing required arguments(s)!")
    missing_arg.each{ log.error("Missing ${it.key}") }
    printhelp = true
}
if (printhelp){
    print(toprint)
    System.exit(0)
}

process publish_pepolar{

    /*
    Map a subject SBRef scan to the derived scan for PEPOLAR
    */

    publishDir path: "${params.out}/${sub}",\
               mode: 'copy',\
               saveAs: { "${name}" },
               overwrite: params.rewrite

    input:
    tuple val(sub), path(img), val(name)

    output:
    path img

    shell:
    '''
    echo "Transferring !{img} to !{name}"
    '''
}

process gen_pepolar{

    /*
    Extract first volume from BOLD to derived scan for PEPOLAR
    */

    label 'FSL'

    input:
    tuple val(sub), val(ser), path(bold)

    output:
    tuple val(sub), val(ser), path("${sub}_${ser}_boldref.nii.gz"),\
    emit: boldref

    shell:
    '''
    #!/bin/bash
    fslroi !{bold} !{sub}_!{ser}_boldref.nii.gz 0 1
    '''
}

// Get full list of available subjects to process
archive_path = "/archive/data/${params.study}/data/nii/"
input_channel = Channel.fromPath("${archive_path}/*", type: 'dir')
                        .map{ i -> i.getBaseName() }

// Only process data in subjects list if provided
if (params.subjects){
    subjects_channel = Channel.fromPath(params.subjects)
                            .splitText() {it.strip()}
    input_channel = input_channel.join(subjects_channel)
}

// Remove subjects that already have an output directory available
if (!params.rewrite){

    out_channel = Channel.fromPath("${params.out}/*", type: 'dir')
                        .map{ o -> [o.getBaseName(), 'o']}
                        .ifEmpty(['', 'o'])
    input_channel = input_channel.join(out_channel, remainder:true)
                                .filter{it.last() == null}
                                .map{i, n -> i}
}

// Pull all BOLD scans
bold_scans = input_channel
                .map{i->[
                            i,
                            params.bold_tags.collect{t ->
                                file("${archive_path}/${i}/*${t}*.nii.gz")
                            }.flatten()
                        ]
                    }
                    .transpose()


// Break scan file into keys
scan_re = /(?<study>[^_]+)_(?<site>[^_]+)_(?<subject>[^_]+)(?<!PHA)_(?<timepoint>[^_]+)_(?<session>[^_]+)/
series_re = /(?<tag>[^_]+)_(?<series>\d+)_(?<description>.*?)(?<ext>.nii.gz)/

filename_re = /${scan_re}_${series_re}/

// Extract series + description + phase encode
pe_re = /(?<=-)(PA|AP)/
bold_sbref = bold_scans
                .map{i,s->[
                    i,s,
                    {
                        def match = s.getName() =~ filename_re
                        match.matches()
                        match
                    }()
                ]}
                .map{i,s,m->[
                    i,s,
                    m.group("series").toInteger(),
                    m.group("description")
                ]}
                .map{i,s,ser,d->[
                    i,s,ser,d,
                    {
                        def match = d =~ pe_re
                        match.find()
                        match[0][0]
                    }(),
                    file("${archive_path}/${i}/${i}*${params.sbref_tag}*"
                    + "${(ser-1).toString().padLeft(2,'0')}*${d}*.nii.gz")
                ]}



workflow{

    main:
    // Generate PEPOLAR from dummy BOLD
    i_gen_pepolar = bold_sbref.filter{it[-1].isEmpty()}
                        .map{i,s,ser,d,pe,f ->[i,ser,s]}
    gen_pepolar(i_gen_pepolar)
    bold2map = bold_sbref.filter{it[-1].isEmpty()}
                    .map{i,s,ser,d,pe,f ->[
                        i,ser,
                        "${i}_${ser.toString().padLeft(2,'0')}_BOLD2FMAP-${pe}.nii.gz"
                    ]}
                    .join(gen_pepolar.out.boldref, by: [0,1])
                    .map{i,ser,n,ref->[
                        i,ref,n
                    ]}


    sbref2map = bold_sbref.filter{!(it[-1].isEmpty())}
                    .map{i,s,ser,d,pe,f ->[
                        i,f,
                        "${i}_${(ser-1).toString().padLeft(2,'0')}"
                        + "_SBREF2FMAP-${pe}.nii.gz"
                    ]}

    publish_pepolar(sbref2map.mix(bold2map))
}
