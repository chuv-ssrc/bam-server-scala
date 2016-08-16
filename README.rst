
A simple BAM server in Scala
============================

This small Play application has been built with the purpose of serving BAM files
for viewing them in IGV.js, while keeping the files secure (since human
genomic data is very sensible).

It relies on a database with a single table and the following schema::

    CREATE TABLE `bam` (
        `id` INT NOT NULL AUTO_INCREMENT PRIMARY KEY,
        `filename` VARCHAR(255) NOT NULL,
        `key` VARCHAR(255) NOT NULL,
        `hash` VARCHAR(255),
        `description` VARCHAR(500)
    ) ENGINE=InnoDB;

In other words, to each BAM file `filename` corresponds a secret `key`.

It is not doing any user authentication. It has no idea which sample each bam file corresponds to.
Instead, it expects the client to know for each sample to read the corresponding secret `key`.
Only this `key` is passed in REST queries.


Configuration
=============

Settings can be edited in `conf/application.conf`. In particular,

- You must indicate ``env.BAM_PATH`` to tell where it can find BAM files for reading.
- It is set to use MySQL. You can change the database driver in the ``db`` section.
- You will probably want to configure CORS in ``play.filters.cors``.

Development
===========

Run the local dev server::

    activator run

Run tests::

    activator test

Build the source (to `target/universal/`)::

    activator dist

Deployment: copy the built archive to destination, decompress, then::

    ./bin/bam-server -v \
        -Dconfig.file=$SETTINGS \
        -Dhttp.port=$PORT \
        -Dhttps.port=$PORT_HTTPS

See also `/scripts/deploy.sh`.

Test that it works::

    ./scripts/test_functional.sh


REST API
========

Let `$key` be the BAM secret identifier in the database.

- GET /

  Says hello, just to test if the server is listening.

- GET /download/index/$key

  Returns the content of the .bai for that BAM file, assuming that the index
  is called <bam_filename>.bai .

- GET /downloadRange/$key

  Returns the content of the BAM file. It is expecting a Range HTTP header
  to extract only the portion of interest, based on the index.
  If suffixed with '.bai', it will return the index just as above.

- GET /download/$key/$region

  Uses samtools (if found) to extract the region (``samtools view -hb <bam> $region``),
  copy it temporarily to ``TEMP_BAM_DIR``, and serve that sub-file.

- GET /read/$key[?region=$region]

  Returns the output of ``samtools view -hb <bam> $region | base64``
  (i.e. need to decode base64 from the response, but at least the binary text can be
  passed using HTTP without creating any temporary file).

- GET /view/$key[?region=$region]

  Returns the reads for the given region in JSON format.
  Currently, it looks like this::

    [
      {
         "chrom": "chr1",
         "name": "long-read-name",
         "start": 1234,
         "end": 45678,
         "cigar": "long-cigar-string"
      },
    ...
    ]

  It can be edited to return whatever you want for you own BAM viewer.
  Look at `htsjdk docs under SAMRecord <https://samtools.github.io/htsjdk/javadoc/htsjdk/>`_
  for available fields.