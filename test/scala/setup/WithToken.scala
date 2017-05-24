package scala.setup

import java.security.PublicKey

import auth.RSA.readPublicKeyFromCertificate
import play.api.test.FakeRequest
import play.api.test.Helpers.GET
import pdi.jwt.{JwtAlgorithm, JwtJson}
import play.api.libs.json.Json



/**
  * Simulate queries from customers with different access rights.
  */
trait WithToken {

  // Public RSA key for test domain jdelafon.eu.auth0.com - bam-server - client ID: G4bFWnpaulRlHScyzznxNZePjNxZG28y
  val auth0PublicKey: PublicKey = readPublicKeyFromCertificate("resources/rsa_keys/auth0.cer")

  // JWTs
  // name = test@test.com (password = 'test'), can access 'sample1' with auth0PublicKey
  val auth0Token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJSUzI1NiIsImtpZCI6Ik1UQTFNRGhHTUVRMVJUUXdOekJGTmpJME9FUkdPRFExUVRneE1rUTBNa0l5TkVWRU56QkROQSJ9.eyJuYW1lIjoidGVzdEB0ZXN0LmNvbSIsImVtYWlsIjoidGVzdEB0ZXN0LmNvbSIsImVtYWlsX3ZlcmlmaWVkIjpmYWxzZSwibmlja25hbWUiOiJ0ZXN0IiwidXNlcl9pZCI6ImF1dGgwfDU5MDlkYWUxNWVhMWM0NWExZWY0ZTBmZCIsImlzcyI6Imh0dHBzOi8vamRlbGFmb24uZXUuYXV0aDAuY29tLyIsInN1YiI6ImF1dGgwfDU5MDlkYWUxNWVhMWM0NWExZWY0ZTBmZCIsImF1ZCI6Ikc0YkZXbnBhdWxSbEhTY3l6em54TlplUGpOeFpHMjh5IiwiZXhwIjozMTQ5NTAzMDQxMCwiaWF0IjoxNDk1MDMwNDEwfQ.cE_Lqum2KNAFKgHMFSrsR-pLrxzoeITVxlmU_xYgPlNo2WzDv9DgCMdcW2V5xI19y0OEo8r8ZLD8YQiR1BLi8bmjyzFTDsxR7dLqSl_Tuv1LWZj4luwk-FkY4iUKUFBVN9MphKJ4cUEDpclPap_VvoMTuF_NYkotl669yYkDxq4F37tN1a9FTf-FB6iE1G5GtGZ4WBVhvcIxAfXGr6lfeW-NytSNsu0XZCUlLWl2MpR1xRWeNbO1IlBFjqlzXB8sHL6cRBHoVgMh8kC4IK6cL9AUHKmaR_8kQBdFBljDnijzDnqP_adDX5VBOS3Cujv16m0BYNokjaMw_C9sUKPNew"
  // name = admin@test.com (password = 'test'), can access 'sample1' and 'sample2' with auth0PublicKey
  val auth0AdminToken = "eyJ0eXAiOiJKV1QiLCJhbGciOiJSUzI1NiIsImtpZCI6Ik1UQTFNRGhHTUVRMVJUUXdOekJGTmpJME9FUkdPRFExUVRneE1rUTBNa0l5TkVWRU56QkROQSJ9.eyJuYW1lIjoiYWRtaW5AdGVzdC5jb20iLCJlbWFpbCI6ImFkbWluQHRlc3QuY29tIiwiZW1haWxfdmVyaWZpZWQiOmZhbHNlLCJuaWNrbmFtZSI6ImFkbWluIiwidXNlcl9pZCI6ImF1dGgwfDU5MWM1OTVkMDNjMDYyMDQyYjVlZDVkMCIsImlzcyI6Imh0dHBzOi8vamRlbGFmb24uZXUuYXV0aDAuY29tLyIsInN1YiI6ImF1dGgwfDU5MWM1OTVkMDNjMDYyMDQyYjVlZDVkMCIsImF1ZCI6Ikc0YkZXbnBhdWxSbEhTY3l6em54TlplUGpOeFpHMjh5IiwiZXhwIjozMTQ5NTAzMDUxMCwiaWF0IjoxNDk1MDMwNTEwfQ.pcu6AdxJKZh7V-YJHQaU_yPjKrrWmp_mOnKjuenU-wXJCjPqrpsXWeLg-YEFp_Trpd1TIb2DvbfOzxubzSpePqn4sMlojwTG8R50e1GGGG7M-omPWZtUCrc-RwHdK2OUp204qs0yRMJnKNkuZsq0Shpis3CwMLfdG1p8xo8NvavgZI5SiS7qeHwM7rD1Y7tjsqZ27DBBovkyiEUVnFnG7C0WS8IuBfEIVxLayXFIakqNcG5b4437mAmeFgt1Jj1m7XVaVuq4k8VW0JVrc0hWt-ujAeFLnbLk7s5Qgx7ky9z7xb2HdmwucTsa2fyhw30aOJirBrs234VVvw9r4BHEqA"
  // name = test@test.com (password = 'test'), from app 2 sharing the same secret:
  private val hmacClaim = Json.obj("iss" -> "testapp2", "name" -> "test@test.com")
  val hmacToken = JwtJson.encode(hmacClaim, "secretHMACkey", JwtAlgorithm.HS256)

  // Respective HTTP headers
  val auth0Header = ("Authorization", "Bearer " + auth0Token)
  val auth0AdminHeader = ("Authorization", "Bearer " + auth0AdminToken)

  def FakeAuthorizedRequest(method: String = GET, url: String) = {
    FakeRequest(method, url).withHeaders(auth0Header)
  }

  def FakeAdminRequest(method: String = GET, url: String) = {
    FakeRequest(method, url).withHeaders(auth0AdminHeader)
  }

  val testUser = "test@test.com"
  val adminUser = "admin@test.com"

  val testSample1 = "sample1"  // samtools view -hb 141.bam chr1:762000-763000
  val testSample2 = "sample2"  // samtools view -hb 141F.recal.bam 1:762000-763000
  val testSample3 = "sample3"  // no bam
  val testSample4 = "sample4"  // no bam
  val notherekey = "notherekey"  // exists in database but not found on disk
  val inactivekey = "inactivekey"  // exists in database but isActive = false

}
