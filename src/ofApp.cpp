#include "ofApp.h"


inline double function_A(double x) {
	return std::exp(-x) * std::sin(8.0 * std::pow(x, 2.0 / 3.0)) + 1.0;
}

inline double function_B(double x) {
	return (14.0 * x - 11.0 * x * x) * std::exp(-2.0 * x);
}
inline double function_B_integral_analytic(double a, double b) {
	auto f = [](double x) {
		return 0.25 * std::exp(-2.0 * x) * (22.0 * x * x - 6.0 * x - 3.0);
	};
	return f(b) - f(a);
}

// Simpson's Rule with cache
template <class Real>
struct SimpsonRange {
	SimpsonRange() {}
	SimpsonRange(std::function<Real(Real)> f, Real a, Real b) :
		_f(f),
		_a(a),
		_fa(f(a)),
		_m((a + b) * Real(0.5)),
		_fm(f(_m)),
		_b(b),
		_fb(f(b))
	{

	}
	SimpsonRange(std::function<Real(Real)> f, Real a, Real fa, Real b, Real fb) :
		_f(f),
		_a(a),
		_fa(fa),
		_m((a + b) * Real(0.5)),
		_fm(f(_m)),
		_b(b),
		_fb(fb)
	{
		_integral = integrate();
	}
	std::tuple<SimpsonRange, SimpsonRange> separate() const {
		return std::make_tuple(
			SimpsonRange(_f, _a, _fa, _m, _fm),
			SimpsonRange(_f, _m, _fm, _b, _fb)
		);
	}

	Real a() const {
		return _a;
	}
	Real fa() const {
		return _fa;
	}
	Real m() const {
		return _m;
	}
	Real fm() const {
		return _fm;
	}
	Real b() const {
		return _b;
	}
	Real fb() const {
		return _fb;
	}
	Real integral() const {
		return _integral;
	}

	Real evaluateApproximation(Real x) const {
		Real h = std::abs(_a - _b) * Real(0.5);
		Real a = (_fa + _fb - Real(2.0) * _fm) / (Real(2.0) * h * h);
		Real b = (_fb - _fa) / (Real(2.0) * h);
		Real c = _fm;
		auto sqr = [](Real x) { return x * x; };
		return a * sqr(x - _m) + b * (x - _m) + c;
	}
private:
	inline Real integrate() {
		return (_b - _a) * Real(1.0 / 6.0) * (_fa + Real(4.0) * _fm + _fb);
	}
	std::function<Real(Real)> _f;
	Real _a;
	Real _fa;
	Real _m;
	Real _fm;
	Real _b;
	Real _fb;
	Real _integral;
};

template <class Real>
Real adaptive_simpson(const SimpsonRange<Real> &simpsonRange, Real eps) {
	SimpsonRange<Real> L;
	SimpsonRange<Real> R;
	std::tie(L, R) = simpsonRange.separate();
	auto L_R_minus_integrate_a_b = L.integral() + R.integral() - simpsonRange.integral();
	if (std::abs(L_R_minus_integrate_a_b) <= Real(15.0) * eps) {
		return L.integral() + R.integral() + L_R_minus_integrate_a_b * Real(1.0 / 15.0);
	}
	return adaptive_simpson(L, eps) + adaptive_simpson(R, eps);
}

template <class Real>
Real adaptive_simpson_visualize(const SimpsonRange<Real> &simpsonRange, Real eps) {
	SimpsonRange<Real> L;
	SimpsonRange<Real> R;
	std::tie(L, R) = simpsonRange.separate();
	auto L_R_minus_integrate_a_b = L.integral() + R.integral() - simpsonRange.integral();
	if (std::abs(L_R_minus_integrate_a_b) <= Real(15.0) * eps) {
		auto draw = [](Real x, Real value) {
			ofSetColor(255);
			ofDrawLine(
				x,
				0,
				x,
				value
			);
			ofDrawCircle(x, value, 0.01);
		};
		draw(L.a(), L.fa());
		ofDrawCircle(L.m(), L.fm(), 0.01);
		draw(L.b(), L.fb());
		draw(R.a(), R.fa());
		ofDrawCircle(R.m(), R.fm(), 0.01);
		draw(R.b(), R.fb());

		auto drawcurve = [](SimpsonRange<Real> simpson) {
			ofPolyline poly;
			int N = 20;
			for (int i = 0; i < N; ++i) {
				float amt = (float)i / (N - 1);
				float x = ofLerp(simpson.a(), simpson.b(), amt);
				float y = simpson.evaluateApproximation(x);
				poly.addVertex(x, y);
			}
			ofSetColor(255, 0, 0);
			poly.draw();
		};
		drawcurve(L);
		drawcurve(R);

		return L.integral() + R.integral() + L_R_minus_integrate_a_b * Real(1.0 / 15.0);
	}
	return adaptive_simpson_visualize(L, eps) + adaptive_simpson_visualize(R, eps);
}

//--------------------------------------------------------------
void ofApp::setup(){

	_camera.setNearClip(0.1f);
	_camera.setFarClip(100.0f);
	_camera.setDistance(5.0f);
}

//--------------------------------------------------------------
void ofApp::update() {
	
}

//--------------------------------------------------------------
void ofApp::draw(){
	//ofEnableDepthTest();

	ofClear(0);
	_camera.begin();
	ofPushMatrix();
	// ofRotateZ(90.0f);
	ofRotateY(90.0f);
	ofSetColor(128);
	ofDrawGridPlane(1.0f);
	ofPopMatrix();

	//ofPushMatrix();
	//ofDrawAxis(50);
	//ofPopMatrix();

	auto func = function_B;

	ofPolyline poly;
	int N = 1000;
	for (int i = 0; i < N; ++i) {
		float x = ofMap(i, 0, N - 1, 0, 2);
		float y = func(x);
		poly.addVertex(x, y);
	}
	ofSetColor(255);
	poly.draw();

	int sampleCount = 0;
	SimpsonRange<double> r([&](double x) {
		sampleCount++;
		return func(x);
	}, 0.0, 2.0);

	double eps = 1.0 * std::pow(0.1, _N);
	double integral = adaptive_simpson_visualize(r, eps);

	_camera.end();

	//ofDisableDepthTest();
	ofSetColor(255);

	{
		char buf[256];
		sprintf(buf, "eps = %.10f", eps);
		ofDrawBitmapString(buf, 10, 10);
	}
	{
		ofDrawBitmapString("Arrow Key <=  =>", 10, 30);
	}

	{
		ofDrawBitmapString("Integral[(14x - 11x^2) exp(-2x), {x, 0, 2}]", 10, 50);
	}

	{
		char buf[256];
		double analytic = function_B_integral_analytic(0.0, 2.0);
		sprintf(buf, "numerical = %.10f, analytic = %.10f, diff = %.10f", integral, analytic, std::abs(integral - analytic));
		ofDrawBitmapString(buf, 10, 70);
	}
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
	if (key == OF_KEY_RIGHT) {
		_N++;
	}
	if (key == OF_KEY_LEFT) {
		_N--;
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
