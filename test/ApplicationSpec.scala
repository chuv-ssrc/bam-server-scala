import akka.stream.Materializer
import org.scalatestplus.play._
import play.api.libs.json._
import play.api.test._
import play.api.test.Helpers._


class ApplicationSpec extends PlaySpec with OneAppPerTest {

  val KEY = "asdf"

  "Routes" should {

    "send 404 on a bad request" in  {
      route(app, FakeRequest(GET, "/fake")).map(status) mustBe Some(NOT_FOUND)
    }

  }

  "HomeController" should {

    "render the index page" in {
      val response = route(app, FakeRequest(GET, "/")).get
      status(response) mustBe OK
      contentType(response) mustBe Some("text/plain")
      contentAsString(response) must include ("BAM server operational.")
    }

  }

  "BamDownloadController" should {

    "return the whole BAM content if no Range is specified" in {
      implicit lazy val materializer: Materializer = app.materializer
      val request = FakeRequest(GET, s"/downloadRange/$KEY")
      val response = route(app, request).get
      contentAsBytes(response).length mustEqual 21794
    }

    "return the whole index if the key has suffix '.bai' in downloadRange" in {
      implicit lazy val materializer: Materializer = app.materializer
      val request = FakeRequest(GET, s"/downloadRange/$KEY.bai")
      val response = route(app, request).get

      status(response) mustBe OK
      contentType(response) mustBe Some(BINARY)
    }

    "return partial BAM content according to the Range header, and status 206" in {
      implicit lazy val materializer: Materializer = app.materializer
      val request = FakeRequest(GET, s"/downloadRange/$KEY")
                    .withHeaders(RANGE -> "bytes=50-150", CONNECTION -> "keep-alive", ACCEPT -> "*/*")
      val response = route(app, request).get
      status(response) mustBe PARTIAL_CONTENT
      contentType(response) mustBe Some(BINARY)
      contentAsBytes(response).length mustEqual 101
    }

    "return the whole BAM content if Range is large enough" in {
      implicit lazy val materializer: Materializer = app.materializer
      val request = FakeRequest(GET, s"/downloadRange/$KEY")
                    .withHeaders(RANGE -> "bytes=0-21793", CONNECTION -> "keep-alive", ACCEPT -> "*/*")
      val response = route(app, request).get
      contentAsBytes(response).length mustEqual 21794
    }

    "return the whole BAM content if Range is too big" in {
      implicit lazy val materializer: Materializer = app.materializer
      val request = FakeRequest(GET, s"/downloadRange/$KEY")
                    .withHeaders(RANGE -> "bytes=0-100000", CONNECTION -> "keep-alive", ACCEPT -> "*/*")
      val response = route(app, request).get
      contentAsBytes(response).length mustEqual 21794
    }

  }


}
