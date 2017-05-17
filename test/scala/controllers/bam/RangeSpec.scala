package scala.controllers.bam

import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._

import scala.setup.WithToken


/**
  * Test RangeController.
  */
class RangeSpec extends PlaySpec with OneAppPerSuite with WithToken {

  val body: JsValue = Json.parse(s"""{"sample": "$testkey"}""")
  val headers = List(
    (AUTHORIZATION -> s"Bearer $auth0Token"),
    (RANGE -> s"bytes=0-65639")
  )

  "RangeController" should {

    "return everything if no range is given (POST)" in {
      val header = (AUTHORIZATION -> s"Bearer $auth0Token")
      val response = route(app, FakeAuthorizedRequest(POST, "/bam/range").withJsonBody(body).withHeaders(header)).get
      status(response) mustBe OK
      contentType(response).get mustBe BINARY
    }

    "provide a slice of the BAM if everything is right (POST)" in {
      val response = route(app, FakeAuthorizedRequest(POST, "/bam/range").withJsonBody(body).withHeaders(headers:_*)).get
      status(response) mustBe PARTIAL_CONTENT
      contentType(response).get mustBe BINARY
    }

    "provide a slice of the BAM if everything is right (GET)" in {
      val response = route(app, FakeAuthorizedRequest(GET, s"/bam/range/$testkey?range=0-65639")).get
      status(response) mustBe PARTIAL_CONTENT
      contentType(response).get mustBe BINARY
    }

    "provide a slice of the BAM if everything is right, with token in URL (GET)" in {
      val response = route(app, FakeRequest(GET, s"/bam/range/$testkey?token=$auth0Token&range=0-65639")).get
      status(response) mustBe PARTIAL_CONTENT
      contentType(response).get mustBe BINARY
    }

    "complain if none of the ranges overlap the resource (POST)" in {
      val headers = List(
        (AUTHORIZATION -> s"Bearer $auth0Token"),
        (RANGE -> s"bytes=151990-355919")
      )
      val response = route(app, FakeAuthorizedRequest(POST, "/bam/range").withJsonBody(body).withHeaders(headers:_*)).get
      status(response) mustBe REQUESTED_RANGE_NOT_SATISFIABLE
    }

    "complain if none of the ranges overlap the resource (GET)" in {
      val response = route(app, FakeAuthorizedRequest(GET, s"/bam/range/$testkey?range=151990-355919")).get
      status(response) mustBe REQUESTED_RANGE_NOT_SATISFIABLE
    }

  }

}