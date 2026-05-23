include {worker} from '../lib/worker'

workflow downstream {
  take:
    mcp
    input

  main:

    // Determine if there are any custom models specified
    inp = Channel.of( mcp.modules['downstream'] )
        .flatten()
        .map{ it -> String m = "${it.name}-model";
            tuple(it, mcp.workflow.containsKey(m) ?
            file(mcp.workflow[m]) : 'built-in') }
        .combine(input)
        .map{ mod, _2, _3 ->
        tuple( '', mod, _2, _3, "${params.in}/downstream/${mod.name}", '') }
    if (mcmicro.Flow.doirun('downstream', mcp.workflow)) {
        worker( mcp, inp, '*.{csv,h5ad}' )
    }

  emit:
    worker.out.res
}
