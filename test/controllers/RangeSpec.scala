package controllers


import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test.Helpers._
import play.api.test._


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
    }

    "provide a slice of the BAM if everything is right (POST)" in {
      val body: JsValue = Json.parse(s"""{"key": "testkey"}""")
      val response = route(app, FakeRequest(POST, "/bam/range").withJsonBody(body).withHeaders(headers:_*)).get
      status(response) mustBe PARTIAL_CONTENT
    }

    "provide a slice of the BAM if everything is right (GET)" in {
      val response = route(app, FakeRequest(GET, s"/bam/range/testkey/0-65639/$token")).get
      status(response) mustBe PARTIAL_CONTENT
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
      val response = route(app, FakeRequest(POST, s"/bam/range/testkey/0-65639/$token")).get
      status(response) mustBe REQUESTED_RANGE_NOT_SATISFIABLE
    }


//    "fail if the key is not known to the database" in {
//      val header = ("Authorization", s"Bearer $token")
//      val body: JsValue = Json.parse(s"""{"key": "xxx"}""")
//      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
//      status(response) mustBe INTERNAL_SERVER_ERROR
//    }
//
//    "fail if the key exists but the bam file is not found" in {
//      val header = ("Authorization", s"Bearer $token")
//      val body: JsValue = Json.parse(s"""{"key": "notherekey"}""")
//      val response = route(app, FakeRequest(POST, "/bai").withJsonBody(body).withHeaders(header)).get
//      status(response) mustBe INTERNAL_SERVER_ERROR
//    }

  }

}