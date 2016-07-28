import akka.stream.Materializer
import org.scalatestplus.play._
import play.api.test._
import play.api.test.Helpers._

/**
 * Add your spec here.
 * You can mock out a whole application including requests, plugins etc.
 * For more information, consult the wiki.
 */
class ApplicationSpec extends PlaySpec with OneAppPerTest {

  val KEY = "asdf"

  "Routes" should {

    "send 404 on a bad request" in  {
      route(app, FakeRequest(GET, "/fake")).map(status) mustBe Some(NOT_FOUND)
    }

  }

  "HomeController" should {

    "render the index page" in {
      val home = route(app, FakeRequest(GET, "/")).get

      status(home) mustBe OK
      contentType(home) mustBe Some("text/plain")
      contentAsString(home) must include ("BAM server operational.")
    }

  }

  "BamDownloadController" should {

    "return the whole BAM content if no Range is specified" in {
      implicit lazy val materializer: Materializer = app.materializer
      val request = FakeRequest(GET, s"/downloadRange/$KEY")
      val indexResponse = route(app, request).get
      contentAsBytes(indexResponse).length mustEqual 21794
    }

    "return the whole index if the key has suffix '.bai' in downloadRange" in {
      implicit lazy val materializer: Materializer = app.materializer
      val request = FakeRequest(GET, s"/downloadRange/$KEY.bai")
      val indexResponse = route(app, request).get

      status(indexResponse) mustBe OK
      contentType(indexResponse) mustBe Some(BINARY)
    }

    "return partial BAM content according to the Range header, and status 206" in {
      implicit lazy val materializer: Materializer = app.materializer
      val request = FakeRequest(GET, s"/downloadRange/$KEY")
                    .withHeaders(RANGE -> "bytes=50-150", CONNECTION -> "keep-alive", ACCEPT -> "*/*")
      val indexResponse = route(app, request).get
      status(indexResponse) mustBe PARTIAL_CONTENT
      contentType(indexResponse) mustBe Some(BINARY)
      contentAsBytes(indexResponse).length mustEqual 101
    }

    "return the whole BAM content if Range is large enough, and status 200" in {
      implicit lazy val materializer: Materializer = app.materializer
      val request = FakeRequest(GET, s"/downloadRange/$KEY")
                    .withHeaders(RANGE -> "bytes=0-21793", CONNECTION -> "keep-alive", ACCEPT -> "*/*")
      val indexResponse = route(app, request).get
      contentAsBytes(indexResponse).length mustEqual 21794
    }

    "return the whole BAM content if Range is too big" in {
      implicit lazy val materializer: Materializer = app.materializer
      val request = FakeRequest(GET, s"/downloadRange/$KEY")
                    .withHeaders(RANGE -> "bytes=0-100000", CONNECTION -> "keep-alive", ACCEPT -> "*/*")
      val indexResponse = route(app, request).get
      contentAsBytes(indexResponse).length mustEqual 21794
    }

  }

  "CountController" should {

    "return an increasing count" in {
      contentAsString(route(app, FakeRequest(GET, "/count")).get) mustBe "0"
      contentAsString(route(app, FakeRequest(GET, "/count")).get) mustBe "1"
      contentAsString(route(app, FakeRequest(GET, "/count")).get) mustBe "2"
    }

  }

}
