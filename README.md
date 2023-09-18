# Structural variant calling from long read sequencing

Workflow to call SVs using different read aligners and calling methods. 
Currently only [NGMLR](https://github.com/philres/ngmlr) and [Sniffles](https://github.com/fritzsedlazeck/Sniffles) implemented.

## Test locally

Quick test with dummy simulated data from [`test`](test).

```sh
miniwdl run --as-me -i test/test.inputs.json workflow.wdl
```
