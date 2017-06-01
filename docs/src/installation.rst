

Installation guidelines
=======================

Bam-server was written in Scala using the Play framework (2.5.4).
Here we summarize the standard deployment guidelines for a Play application.
For more details, refer to the
`Play docs <https://www.playframework.com/documentation/2.5.x/Production>`_.


Requirements
------------

- Java RE 8+

To work from source:

- Scala 2.11.7+
- sbt 0.13.11+


Running in dev mode
-------------------

Running the development server can be useful either to quickly test the application
before setting up a production installation, or when you want to edit the source code.
From source, in the root directory:

.. code:: bash

    sbt run

By default it gets served at `localhost:9000` (try in your browser's address bar).

To use the test configuration instead of the default (recommended):

.. code:: bash

  	sbt run -Dconfig.resource=application.test.conf

To run the tests:

.. code:: bash

    sbt test

Tests use a built-in H2 database with hard-coded data that gets destroyed and recreated on each run.
Its contents can be found in `/conf/evolutions/default/1.sql`.


Deployment
----------

Get the latest release `from GitHub <https://github.com/jdelafon/bam-server-scala/releases/latest>`_.
A precompiled release for v1.0 is available
`here <https://github.com/jdelafon/bam-server-scala/releases/download/v1.0/bam-server-1.0-SNAPSHOT.zip>`_,
or you can `build from source`_.
Move the build archive to destination, decompress, enter the directory, then launch the server
(`netty <https://netty.io/>`_):

.. code:: bash

    ./bin/bam-server

By default it gets served at `localhost:9000` (try in your browser's address bar).
When it starts, a file `RUNNING_PID` is created in the root directory that contains the ID
of the running process. To stop the application from running, kill this process
and do not forget to remove the file:

.. code:: bash

    cat RUNNING_PID | xargs kill -9
    rm RUNNING_PID

You can configure a proxy to make it available from the outside world. For example with Apache:

.. code:: apache

    <VirtualHost *:443>
        ...
        ProxyPass         /bamserver  http://localhost:9000
        ProxyPassReverse  /bamserver  http://localhost:9000
        ...
    </Virtualhost>

The service becomes available at ``https://<host>/bamserver/``.
Apache can take care itself of protecting external transactions with HTTPS.
To use another proxy than Apache, refer to
`Play HTTPServer docs <https://www.playframework.com/documentation/2.5.x/HTTPServer>`_.


.. _Configuration:

Configuration
-------------

When deploying, you will certainly want to configure at least the server port, the application secret,
or completely replace the configuration file (recommended), for instance:

.. code:: bash

    ./bin/bam-server -v \
        -Dconfig.file=../application.prod.conf \
        -Dhttp.port=8870 \
        -Dhttps.port=8871

Settings can be edited in `/conf/application.conf`, but we recommend specifying a new config file
with ``-Dconfig.file=<path/to/config>`` as shown above.

All you want to know about configuration files and options for Play is described
`there <https://www.playframework.com/documentation/2.5.x/Configuration>`_.

For bam-server in particular,

- You must indicate ``env.BAM_PATH`` to tell where it can find BAM files for reading.
- It is set to use SQLite. You can change the database settings in the ``db`` section.
- You will probably want to configure CORS in ``play.filters.cors`` to white-list your client
  and restrict others.

A default SQLite database is provided with the project source (`/resources/db/bam-server`),
already with the schema and a few test data.


.. _build from source:

Building from source
--------------------

.. code:: bash

    git clone --depth 1 https://github.com/jdelafon/bam-server-scala.git
    cd bam-server-scala/
    sbt dist

It creates a distribution archive under `/target/universal/`.

