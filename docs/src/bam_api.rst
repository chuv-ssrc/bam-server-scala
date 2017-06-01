

BAM query API
=============

Bam-server can return portions of BAM files using 3 different strategies:

1. `Extracting a range of bytes`_ (IGV.js uses that);
2. `Using samtools`_ (if available);
3. `Using Picard-tools/htsjdk`_ to return reads in JSON format.


Authorization
-------------

All requests expect an Authorization header to be added to the request,
containing a signed Bearer token (JWT). For example:

.. code:: bash

    curl -i -H "Authorization: Bearer xxxx.yyyy.zzzz" -X POST -d '{"sample": "SAMPLE1"} http://localhost:9000/bai
    curl -i -H "Authorization: Bearer xxxx.yyyy.zzzz" http://localhost:9000/bai/SAMPLE1

Alternatively, all request can accept a URL argument `?token=<JWT>` to pass the token.
This mode is not recommended but sometimes necessary, depending on client library constraints.
For example:

.. code:: bash

    curl -i http://localhost:9000/bai/SAMPLE1?token=xxxx.yyyy.zzzz


API
---

Let `:sample`/`<sample>` be a sample identifier in the database;

Let `<range>` be a bytes range, formatted as "123-456";

Let `<region>` be a genomic region, formatted as "chr1:10000-20000".
The chromosome names are those referenced in the BAM file header.


.. code:: bash

    GET /

Says "BAM server operational.", just to test if the server is listening.


BAM index
.........

.. code:: bash

    POST /bai                                   # data = '{"sample": "<sample>"}'
    GET /bai/:sample

Returns the content of the index (.bai) for that BAM file.


.. _Extracting a range of bytes:

Range query
...........

.. code:: bash

    POST /bam/range                             # data = '{"sample": "<sample>"}'
    GET /bam/range/:sample

Returns the content of the BAM file, expecting an
HTTP `Range <https://developer.mozilla.org/en-US/docs/Web/HTTP/Range_requests>`_ header
to extract only the bytes range of interest - likely based on the BAM index.
The bytes range can also be passed as an argument to the request (`?range=<range>`).
Return the content as binary.


.. _Using samtools:

Using samtools
..............

.. code:: bash

    POST /bam/samtools?region=<region>          # data = '{"sample": "<sample>"}'
    GET /bam/samtools/:sample?region=<region>

Uses samtools (if available) to extract the region (``samtools view -hb <bam> <region>``).
Return the content as binary.


.. _Using Picard-tools/htsjdk:

Reads in JSON format
....................

.. code:: bash

    POST /bam/json?region=<region>              # data = '{"sample": "<sample>"}'
    GET /bam/json/:sample?region=<region>

Returns the reads for the given region in JSON format,
using the `htsjdk <http://samtools.github.io/htsjdk/>`_ library.
The fields correspond to the SAM file columns:

.. code:: bash

    [
      {
         "name": "HISEQ:206:C8E95ANXX:3:2113:2451:6639",   // read name
         "flag": 99,
         "chrom": "chr1",       // reference name
         "start": 1234,         // leftmost mapping position
         "end": 1334,           // rightmost mapping position
         "mapq": 50,            // mapping quality
         "cigar": "101M",       // cigar string
         "rnext": "=",          //
         "pnext": 4567,         //
         "tlen": 283,           // template length, aka insert size
         "seq": "AATTAGGA...",  // [ACGTN=.],
         "qual": "AB<B@G>F..."  // per-base quality
      },
    ...
    ]

