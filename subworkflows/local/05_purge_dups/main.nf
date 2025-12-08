/*
 * Workflow based around the purge_dups tool
 * https://github.com/dfguan/purge_dups
 */

include { combineByMetaKeys        } from "../../../modules/local/functions"
include { constructAssemblyRecord  } from "../../../modules/local/functions"
include { getEachAssembly          } from "../../../modules/local/functions"
include { ITERATIVE_PURGEDUPS      } from "../../../modules/local/purgedups/main"

workflow PURGE_DUPLICATES {

    take:
    ch_assemblies // [ meta, assembly ]
    ch_hifi       // [ meta, hifi_merged ]
    ch_kmer_cov   // [ meta, kmer_cov ]

    main:
    // Update with k-mer coverage
    ch_updated_hifi = combineByMetaKeys(
        ch_hifi,
        ch_kmer_cov,
        keySet: ['id','sample'],
        meta: 'lhs'
    ).map { meta, reads, kmer_cov ->
        tuple( kmer_cov ? meta + [kmercov: kmer_cov] : meta, reads )
    }
    // Combine reads with assemblies and prepare for alignment
    ch_reads_with_assemblies = combineByMetaKeys(
        ch_updated_hifi,
        getEachAssembly(ch_assemblies),
        keySet: ['id', 'sample'],
        meta: 'rhs'
    )
    ITERATIVE_PURGEDUPS( ch_reads_with_assemblies )

    ch_purged_assemblies = constructAssemblyRecord(
        ITERATIVE_PURGEDUPS.out.purged.transpose(),
        true
    )

    logs = ITERATIVE_PURGEDUPS.out.logs
        .flatMap { _meta, file -> file }

    versions = ITERATIVE_PURGEDUPS.out.versions

    emit:
    assemblies = ch_purged_assemblies
    logs
    versions
}