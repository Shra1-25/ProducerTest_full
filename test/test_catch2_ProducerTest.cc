/*#include "catch.hpp"
#include "\\home\\cmsusr\\CMSSW_10_6_8\\src\\demo\\plugins\\ProducerTest.h"
#include "FWCore/TestProcessor/interface/TestProcessor.h"
#include "FWCore/Utilities/interface/Exception.h"

int test_tf(){
  ProducerTest pred=new ProducerTest();
  pred.produce();
  return 1;
}

TEST_CASE("Testing sum and product","[sum_product]"){
  REQUIRE(test_tf()==1);
}
static constexpr auto s_tag = "[ProducerTest]";

TEST_CASE("Standard checks of ProducerTest", s_tag) {
  const std::string baseConfig{
R"_(from FWCore.TestProcessor.TestProcess import *
process = TestProcess()
process.toTest = cms.EDProducer("ProducerTest"
#necessary configuration parameters
 )
process.moduleToTest(process.toTest)
)_"
  };
  
  edm::test::TestProcessor::Config config{ baseConfig };  
  SECTION("base configuration is OK") {
    REQUIRE_NOTHROW(edm::test::TestProcessor(config));
  }
  
  SECTION("No event data") {
    edm::test::TestProcessor tester(config);
    
    REQUIRE_THROWS_AS(tester.test(), cms::Exception);
    //If the module does not throw when given no data, substitute 
    //REQUIRE_NOTHROW for REQUIRE_THROWS_AS
  }
  
  SECTION("beginJob and endJob only") {
    edm::test::TestProcessor tester(config);
    
    REQUIRE_NOTHROW(tester.testBeginAndEndJobOnly());
  }

  SECTION("Run with no LuminosityBlocks") {
    edm::test::TestProcessor tester(config);
    
    REQUIRE_NOTHROW(tester.testRunWithNoLuminosityBlocks());
  }

  SECTION("LuminosityBlock with no Events") {
    edm::test::TestProcessor tester(config);
    
    REQUIRE_NOTHROW(tester.testLuminosityBlockWithNoEvents());
  }

}

//Add additional TEST_CASEs to exercise the modules capabilities
*/
