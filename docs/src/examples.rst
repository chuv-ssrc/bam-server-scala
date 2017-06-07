

Implementation examples
=======================

Here are some practical solutions using external products that we have used with success
in practice. We do not advertise for them but only provide these resources to illustrate
concretely how bam-server can be used.


Using IGV.js as reads viewer
----------------------------

The IGV team has produced a Javascript version of its famous BAM viewer, available at
`<https://github.com/igvteam/igv.js/tree/master>`_.

All you need to do to use it is pass the right URLs to its ``options.tracks`` parameter:

.. code:: javascript

    options = {
        showNavigation: true,
        showRuler: true,
        genome: "hg19",
        locus: "chr1:761997-762551"
        tracks: [
            {
                url: 'http://localhost:9000/bam/range/'+ sampleName +'?token='+ AuthService.getToken(),
                indexURL: 'http://localhost:9000/bai/'+ sampleName +'?token='+ AuthService.getToken(),
                type: "alignment",
                name: sampleName + ' [Protected resource!!]',
                displayMode: 'SQUISHED',
                alignmentRowHeight: 4,
            }
        ]
    };

Here we use

.. code:: javascript

    /bai/:sample?token=<token>
    /bam/range/:sample?token=<token>

from the bam-server API to read respectively the BAM index and the BAM itself.
IGV.js is able to parse the index to find which bytes range to ask for in a RANGE request.

.. note::

    The IGV.js team is about to release a way to pass Authorization headers
    instead of having to pass the token in the URL (which can be unsafe with reusable tokens).


Result:

.. figure:: /images/demo.png
   :width: 100%
   :alt: Demo interface using bam-server to display read alignments

   Demo interface using bam-server to display read alignments.


Example of using Auth0 as authorization server
----------------------------------------------

`Auth0 <https://manage.auth0.com>`_ is a commercial solution that provides an easy-to-manage authorization server
(for free if you have few users and minimal requirements) implementing the
`OAuth2 <https://oauth.net/2/>`_ protocol.
Its only role is to identify users of the client app and return signed tokens.


1. Create an Auth0 client and register it to bam-server

    1. Create an Auth0 account

    2. Create a new client app, say "bam-server". It is given a Client ID and a Client secret.

    3. Go to "Clients > Show advanced settings". In the "OAuth" tab, choose "RS256". This changes
       the encryption protocol from the default HMAC (symmetric, shared secret)
       to RSA (asymmetric private/public - more secure).

    4. Still in the advanced settings, copy the public certificate from "Certificates > Signing Certificate"
       in CER or PEM format.

    5. Create an app in the `apps` table of the bam-server database, with the "Domain name" in the `iss` column.

    6. Replace line returns in the public certificate by ``\n`` and paste this as a single line in the `key` column.


2. Obtain a token (for testing; the client app is out of this scope)

    1. Let's say your client web app is served as localhost:3000.

    2. Add "localhost:3000" to Auth0 "Clients > Allowed callback URLs".

    3. Point your browser to::

           https://<domain>/authorize
               ?scope=openid name nickname email user_id
               &response_type=token
               &client_id=<client_id>
               &redirect_uri=localhost:3000
               &nonce=1234&state=1234

       where you replace `<domain>` and `<client_id>` by your own Auth0 client settings.
       N.B. The "scope" argument is what permits to have the fields "name", "email", etc. in the token,
       which we can use a user identifier.

    4. Log in your client app; you get redirected to localhost:3000/ and you can find an "id_token" in the url.
       This is a JWT that bam-server accepts (*not* "access_token"). For instance::

           http://localhost:3000?access_token=ZVvVzzBIucYdhnW3&expires_in=86400&id_token=xxxxxxx.yyyyyyy.zzzzzzz&token_type=Bearer&state=1234

       To decode the body of the token (which is public!), you can use `jwt.io <https://jwt.io/>`_.
       Note the `iss` and `name` claims.


3. Create a user, a sample, a users_samples entry in the database

    1. Create a new user using the `name` from the JWT as `username`,
       and the ID of the app you created in 1. as `app_id`.

       .. code:: SQL

           INSERT INTO users(app_id, username) VALUES (<app_id>, <username>);

    2. Add a sample:

       .. code:: SQL

           INSERT INTO bam(sample, filename) VALUES (<sample_name>, <bam_filename>);

       Make sure the BAM file is available at `env.BAM_PATH` (see Configuration in :doc:`./installation`).

    3. Attribute the sample to a user:

       .. code:: SQL

           INSERT INTO users_bam(user_id, bam_id) VALUES (<user_id>, <bam_id>);

       where `<user_id>` and `<bam_id>` are the respective ids of the rows inserted above.


4. Use the client app to query

   Try to read a BAM index on behalf of the user you created
   (i.e. using the token we got in 2.):

       curl -i -H "Authorization: Bearer xxxxxxx.yyyyyyy.zzzzzzz" http://localhost:9000/bai/SAMPLE1



