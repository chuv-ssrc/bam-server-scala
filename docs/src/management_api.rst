

Management API
==============

A mini REST API has also been added to facilitate the entry of apps, users, samples
and users_samples attributions in the database.

Only registered users tagged as `isAdmin` in the database have access to these functions.
Like BAM query routes, these expect an Authorization header with a signed Bearer token
to identify the requester (see :doc:`./auth`).


Add/remove an app
.................

To register an application to the bam-server,
i.e. make bam-server understand (verify) signed tokens coming from the registered applications.

.. code:: bash

    PUT /apps

Expects a JSON body of the type

.. code:: json

    {
      "iss": "my-app",
      "key": "secret",
      "description": ""
    }

- **iss**: the issuer, the app identifier such as present in JWTs under the `iss` claim.
- **key**: the signature verification key. Either a shared HMAC secret or an RSA public certificate
  (.pem or .cer formats, replacing newlines by ``\n``).
- **description**: [optional] a description of the app.

.. code:: bash

    DELETE /apps

Expects a JSON body of the type

.. code:: json

    {
      "iss": "my-app",
    }

.. note::

    One cannot add an app that already exists (same **iss**).


Add/remove users
................

To register users that can make BAM queries.
Users will be given the same `appId` as the admin inserting them.

.. code:: bash

    PUT /users

Expects a JSON body of the type

.. code:: json

    {
      "users": [
        {"username": "user1"},
        {"username": "user2"}
      ]
    }

.. code:: bash

    DELETE /users

Expects a JSON body of the type

.. code:: json

    {
      "users": ["user1", "user2"],
    }

.. note::

    One cannot add users that already exist, or delete users that do not exist.
    If that happens for one of the users of the query, nothing is inserted or deleted at all.


Add/remove samples
..................

To register BAM files available for query.

.. code:: bash

    PUT /samples

Expects a JSON body of the type

.. code:: json

    {
      "samples": [
        {
          "name": "A",
          "filename": "/"
        },
        {
          "name": "B",
          "filename": "/"
        }
      ]
    }

- **name**: the sample identifer, such as used in the BAM query API.
- **filename**: file name, or path to the BAM file relatively to the configured `BAM_PATH`.

.. code:: bash

    DELETE /samples

Expects a JSON body of the type

.. code:: json

    {
      "samples": ["A", "B"]
    }

.. note::

    One cannot add samples that already exist, or delete samples that do not exist.
    If that happens for one of the samples of the query, nothing is inserted or deleted at all.


Add/remove an attribution
.........................

To give or revoke access of a certain user to a certain BAM file.

.. code:: bash

    PUT /users_samples
    DELETE /users_samples

Expect a JSON body of the type

.. code:: json

    {
      "users_samples": [
        {
          "sample": "S1",
          "username": "A"
        },
        {
          "sample": "S2",
          "username": "B"
        }
      ]
    }

.. note::

    All user identifiers and sample identifiers in a query must be found in the database.
    If it is not the case of one of them, nothing is inserted or deleted at all.
