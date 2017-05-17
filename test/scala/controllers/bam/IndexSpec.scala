package scala.controllers.bam

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._

import scala.setup.WithToken


/**
  * Test IndexController.
  */
class IndexSpec extends PlaySpec with OneAppPerSuite with WithToken {

  /*
    This first part should certainly be tested outside of the IndexController scope. For instance:

      import org.scalatest.mock.MockitoSugar
      class BamQuerySpec extends PlaySpec with MockitoSugar {
        val ctrl = new BamQueryController(mock[Database])
      }
   */

  "`parseBamRequestFromPost`" should {

    "fail if there is no JSON body" in {
      val response = route(app, FakeAuthorizedRequest(POST, "/bai")).get
      status(response) mustBe UNSUPPORTED_MEDIA_TYPE
      // Caught at the 'Action(parse.json)' level
    }

    "fail if there is no key 'sample' in the body" in {
      val body: JsValue = Json.parse(s"""{"yyyy": "xxx"}""")
      val response = route(app, FakeAuthorizedRequest(POST, "/bai").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("No key found in request body")
    }

  }

  "`sampleNameToBamRequest`" should {

    "fail if the key is not known to the database (POST)" in {
      val body: JsValue = Json.parse(s"""{"sample": "xxx"}""")
      val response = route(app, FakeAuthorizedRequest(POST, "/bai").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("No corresponding BAM file")
    }

    "fail if the key exists but the bam file is not found on disk (POST)" in {
      val body: JsValue = Json.parse(s"""{"sample": "$notherekey"}""")
      val response = route(app, FakeAuthorizedRequest(POST, "/bai").withJsonBody(body)).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("This BAM file cannot be found")
    }

    "fail if the key is not known to the database (GET)" in {
      val response = route(app, FakeAuthorizedRequest(GET, s"/bai/xxx?token=$auth0Token")).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("No corresponding BAM file")
    }

    "fail if the key exists but the bam file is not found (GET)" in {
      val response = route(app, FakeAuthorizedRequest(GET, s"/bai/$notherekey?token=$auth0Token")).get
      status(response) mustBe INTERNAL_SERVER_ERROR
      contentAsString(response) must startWith("This BAM file cannot be found")
    }

  }

  "IndexController" should {

    "provide the BAM index if everything is right (POST)" in {
      val body: JsValue = Json.parse(s"""{"sample": "$testkey"}""")
      val response = route(app, FakeAuthorizedRequest(POST, "/bai").withJsonBody(body)).get
      status(response) mustBe OK
      contentType(response).get mustBe BINARY  // application/octet-stream
    }

    "provide the BAM index if everything is right (GET)" in {
      val response = route(app, FakeAuthorizedRequest(GET, s"/bai/$testkey")).get
      status(response) mustBe OK
      contentType(response).get mustBe BINARY
    }

    "provide the BAM index if everything is right, with token in URL (GET)" in {
      val response = route(app, FakeRequest(GET, s"/bai/$testkey?token=$auth0Token")).get
      status(response) mustBe OK
      contentType(response).get mustBe BINARY
    }

  }

}
