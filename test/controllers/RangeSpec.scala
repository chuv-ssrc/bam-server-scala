package controllers


import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._


/**
  * Test RangeController.
  */
class RangeSpec extends PlaySpec with OneAppPerSuite {

  val token = "asdf"
  val headers = List(
    (AUTHORIZATION -> s"Bearer $token"),
    (RANGE -> s"bytes=0-65639")
  )

  "RangeController" should {

    "Return everything if no range is given (POST)" in {
      val body: JsValue = Json.parse(s"""{"key": "testkey"}""")
      val headers = (AUTHORIZATION -> s"Bearer $token")
      val response = route(app, FakeRequest(POST, "/bam/range").withJsonBody(body).withHeaders(headers)).get
      status(response) mustBe OK
      contentType(response).get mustBe BINARY
    }

    "provide a slice of the BAM if everything is right (POST)" in {
      val body: JsValue = Json.parse(s"""{"key": "testkey"}""")
      val response = route(app, FakeRequest(POST, "/bam/range").withJsonBody(body).withHeaders(headers:_*)).get
      status(response) mustBe PARTIAL_CONTENT
      contentType(response).get mustBe BINARY
    }

    "provide a slice of the BAM if everything is right (GET)" in {
      val response = route(app, FakeRequest(GET, s"/bam/range/testkey/$token?range=0-65639")).get
      status(response) mustBe PARTIAL_CONTENT
      contentType(response).get mustBe BINARY
    }

    "Complain if none of the ranges overlap the resource (POST)" in {
      val body: JsValue = Json.parse(s"""{"key": "testkey"}""")
      val headers = List(
        (AUTHORIZATION -> s"Bearer $token"),
        (RANGE -> s"bytes=151990-355919")
      )
      val response = route(app, FakeRequest(POST, "/bam/range").withJsonBody(body).withHeaders(headers:_*)).get
      status(response) mustBe REQUESTED_RANGE_NOT_SATISFIABLE
    }

    "Complain if none of the ranges overlap the resource (GET)" in {
      val response = route(app, FakeRequest(GET, s"/bam/range/testkey/$token?range=151990-355919")).get
      status(response) mustBe REQUESTED_RANGE_NOT_SATISFIABLE
    }

  }

}