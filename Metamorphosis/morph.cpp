/****************************************************************************
 * THE GRAND METAMORPHOSIS
 * CS148 Assignment #4 - Fall 2010, Stanford University
 
 Student: Daniel Capo
 ****************************************************************************/

#include "st.h"
#include "stglut.h"
#include "parseConfig.h"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

// --------------------------------------------------------------------------
// Structure to contain an image feature for a morph. A feature is a directed
// line segment from P to Q, with coordinates in pixel units relative to the
// lower-left corner of the image.
// --------------------------------------------------------------------------

struct Feature
{
    STPoint2 P, Q;
    Feature(const STPoint2 &p, const STPoint2 &q) : P(p), Q(q) { }
};

// --------------------------------------------------------------------------
// Constants, a few global variables, and function prototypes
// --------------------------------------------------------------------------

const int kWindowWidth  = 512;
const int kWindowHeight = 512;
const int kFrames       = 30;   // number of frames to generate

STImage *gDisplayedImage = 0;   // an image to display (for testing/debugging)

std::vector<Feature> gSourceFeatures;   // feature set on source image
std::vector<Feature> gTargetFeatures;   // corresponding features on target

// Copies an image into the global image for display
void DisplayImage(STImage *image);

// --------------------------------------------------------------------------
// CS148 TODO: Implement the functions below to compute the morph
// --------------------------------------------------------------------------


float Lerp(float c1, float c2, float t) {
    return c1 + t * (c2 - c1);
}

STColor4ub colorLerp(STColor4ub c1, STColor4ub c2, float t) {
    float r = Lerp(c1.r, c2.r, t);
    float g = Lerp(c1.g, c2.g, t);
    float b = Lerp(c1.b, c2.b, t);
    STColor4ub result(r,g,b);
    return result;
}

STColor4ub biLerp(STPoint2& X_prime, STImage *image) {
    
    STPoint2 v0(floorf(X_prime.x), floorf(X_prime.y));
    STPoint2 v1(ceilf(X_prime.x), floorf(X_prime.y));
    STPoint2 v2(floorf(X_prime.x), ceilf(X_prime.y));
    STPoint2 v3(ceilf(X_prime.x), ceilf(X_prime.y));
    STColor4ub v0C, v1C, v2C, v3C; 
    
    // Edge cases
    if (v0.x >= image->GetWidth() || v0.y >= image->GetHeight())  v0C = STColor4ub(0,0,0);
    else v0C = STColor4ub(image->GetPixel(v0.x, v0.y));
    if (v1.x >= image->GetWidth() || v1.y >= image->GetHeight())  v1C = STColor4ub(0,0,0);
    else v1C = STColor4ub(image->GetPixel(v1.x, v1.y));
    if (v2.x >= image->GetWidth() || v2.y >= image->GetHeight())  v2C = STColor4ub(0,0,0);
    else v2C = STColor4ub(image->GetPixel(v2.x, v2.y)); 
    if (v3.x >= image->GetWidth() || v3.y >= image->GetHeight())  v3C = STColor4ub(0,0,0);
    else v3C = STColor4ub(image->GetPixel(v3.x, v3.y));
    
    float s = X_prime.x - v0.x;
    float t = X_prime.y - v0.y;
    STColor4ub v01C(colorLerp(v0C, v1C, s));
    STColor4ub v23C(colorLerp(v2C, v3C, s));
    STColor4ub v(colorLerp(v01C, v23C, t));
    return v;
}

/**
 * Compute a linear blend of the pixel colors in two provided images according
 * to a parameter t.
 */
STImage *BlendImages(STImage *image1, STImage *image2, float t)
{
    int minWidth = std::min(image1->GetWidth(), image2->GetWidth());
    int minHeight = std::min(image1->GetHeight(), image2->GetHeight());
    STImage *result = new STImage(minWidth, minHeight);
    for (int x = 0; x < minWidth; x++) {
        for (int y = 0; y < minHeight; y++) {
            STColor4ub p1 = image1->GetPixel(x, y);
            STColor4ub p2 = image2->GetPixel(x, y);
            STColor4ub resultPixel(colorLerp(p1, p2, t));
            result->SetPixel(x, y, resultPixel);
        }
    }
    return result;
}


/**
 * Compute a field morph on an image using two sets of corresponding features
 * according to a parameter t.  Arguments a, b, and p are weighting parameters
 * for the field morph, as described in Beier & Nelly 1992, section 3.
 */
STImage *FieldMorph(STImage *image,
                    const std::vector<Feature> &sourceFeatures,
                    const std::vector<Feature> &targetFeatures,
                    float t, float a, float b, float p)
{
    STImage *result = new STImage(image->GetWidth(), image->GetHeight());
    for (int x = 0; x < result->GetWidth(); x++) {
        for (int y = 0; y < result->GetHeight(); y++) {
            STPoint2 X(x,y);
            STVector2 dSum(0,0);
            float weightSum = 0;
            for (int i = 0; i < targetFeatures.size(); i++) {
                float Pi_x = Lerp(sourceFeatures[i].P.x, targetFeatures[i].P.x, t);
                float Pi_y = Lerp(sourceFeatures[i].P.y, targetFeatures[i].P.y, t);
                float Qi_x = Lerp(sourceFeatures[i].Q.x, targetFeatures[i].Q.x, t); 
                float Qi_y = Lerp(sourceFeatures[i].Q.y, targetFeatures[i].Q.y, t); 

                STPoint2 Pi(Pi_x, Pi_y);
                STPoint2 Qi(Qi_x, Qi_y);
                STVector2 PiQi(Qi - Pi);
                STVector2 PiX(X - Pi);
                STVector2 perpPiQi(-PiQi.y, PiQi.x);
                float u = STVector2::Dot(PiX, PiQi) / PiQi.LengthSq();
                float v = STVector2::Dot(PiX, perpPiQi) / PiQi.Length();
                
                STPoint2 Pi_prime = sourceFeatures[i].P;
                STPoint2 Qi_prime = sourceFeatures[i].Q;
                STVector2 PiQi_prime(Qi_prime - Pi_prime);
                STVector2 perpPiQi_prime(-PiQi_prime.y, PiQi_prime.x);
                STPoint2 Xi_prime( Pi_prime + (u * PiQi_prime) + (v * perpPiQi_prime)/PiQi_prime.Length() );
                
                STVector2 Di(Xi_prime - X);
                float dist;
                if (u < 0) dist = STPoint2::Dist(Pi, X);
                else if (u > 1) dist = STPoint2::Dist(Qi,X);
                else dist = abs(v);
                
                float weight = powf( (powf(PiQi.Length(),p) / (a + dist)), b);
                dSum += Di * weight;
                weightSum += weight;
            }
            STPoint2 X_prime = X + dSum / weightSum;
            if (   X_prime.x < 0 
                || X_prime.x >= image->GetWidth()
                || X_prime.y < 0 
                || X_prime.y >= image->GetHeight()) {
                // set pixel to be WHITE
               // result->SetPixel(x, y, STColor4ub(255,255,255,255));
            } else {
                STColor4ub X_prime_color = biLerp(X_prime, image);
                result->SetPixel(x, y, X_prime_color);
            }
        }
    }
    return result;
}

/**
 * Compute a morph between two images by first distorting each toward the
 * other, then combining the results with a blend operation.
 */
STImage *MorphImages(STImage *sourceImage, const std::vector<Feature> &sourceFeatures,
                     STImage *targetImage, const std::vector<Feature> &targetFeatures,
                     float t, float a, float b, float p)
{
    STImage *image1 = FieldMorph(sourceImage, sourceFeatures, targetFeatures, t, a, b, p);
    STImage *image2 = FieldMorph(targetImage, targetFeatures, sourceFeatures, 1-t, a, b, p);
    
    STImage *result = BlendImages(image1, image2, t);
    return result;
}

/**
 * Compute a morph through time by generating appropriate values of t and
 * repeatedly calling MorphImages(). Saves the image sequence to disk.
 */
void GenerateMorphFrames(STImage *sourceImage, const std::vector<Feature> &sourceFeatures,
                         STImage *targetImage, const std::vector<Feature> &targetFeatures,
                         float a, float b, float p)
{
    // iterate and generate each required frame
    float t = 0;
    for (int i = 0; i <= kFrames; ++i)
    {
        std::cout << "Metamorphosizing frame #" << i << "...";
        float ease_t = powf(t, 2.f)*(3-2*t);
        STImage *result = MorphImages(sourceImage, sourceFeatures, targetImage, targetFeatures, ease_t, a, b, p);
        t += (1.0/30.0);
        // generate a file name to save
        std::ostringstream oss;
        oss << "frame" << std::setw(3) << std::setfill('0') << i << ".png";

        // write and deallocate the morphed image
        if (result) {
            result->Save(oss.str());
            delete result;
        }

        std::cout << " done." << std::endl;
    }
}

// --------------------------------------------------------------------------
// Utility and support code below that you do not need to modify
// --------------------------------------------------------------------------

/**
 * Copies an image into the global image for display
 */
void DisplayImage(STImage *image)
{
    // clean up the previous image
    if (gDisplayedImage) {
        delete gDisplayedImage;
        gDisplayedImage = 0;
    }

    // allocate a new image and copy it over
    if (image) {
        gDisplayedImage = new STImage(image->GetWidth(), image->GetHeight());
        size_t bytes = image->GetWidth() * image->GetHeight() * sizeof(STImage::Pixel);
        memcpy(gDisplayedImage->GetPixels(), image->GetPixels(), bytes);
    }
}

/**
 * Display callback function draws a single image to help debug
 */
void DisplayCallback()
{
    glClearColor(.2f, 2.f, 2.f, 1.f);
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if (gDisplayedImage)
        gDisplayedImage->Draw();

    glutSwapBuffers();
}

/**
 * Window resize callback function
 */
void ReshapeCallback(int w, int h)
{
    glViewport(0, 0, w, h);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0, w, 0, h);
}

/**
 * Keyboard callback function
 */
void KeyboardCallback(unsigned char key, int x, int y)
{
    switch (key)
    {
        // exit program on escape press
        case 27:
            exit(0);
            break;
        // save the currently displayed image if S is pressed
        case 's':
        case 'S':
            if (gDisplayedImage)
                gDisplayedImage->Save("screenshot.png");
            break;
        default:
            break;
    }
}

/**
 * This function is called by the parsing functions to populate the feature sets
 */
void AddFeatureCallback(STPoint2 p, STPoint2 q, ImageChoice image)
{
    if (image == IMAGE_1 || image == BOTH_IMAGES)
        gSourceFeatures.push_back(Feature(p, q));
    if (image == IMAGE_2 || image == BOTH_IMAGES)
        gTargetFeatures.push_back(Feature(p, q));
}

/**
 * Program entry point
 */
int main(int argc, char* argv[])
{
    glutInit(&argc, argv);
    glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGB );
    glutInitWindowPosition(20, 20);
    glutInitWindowSize(kWindowWidth, kWindowHeight);
    glutCreateWindow("Metamorphosis: CS148 Assignment 4");

    glutDisplayFunc(DisplayCallback);
    glutReshapeFunc(ReshapeCallback);
    glutKeyboardFunc(KeyboardCallback);

    //
    // load the configuration from config.txt, or other file as specified
    //
    std::string configFile = "config.txt";
    if (argc > 1) configFile = argv[1];

    char sourceName[64], targetName[64];
    char saveName[64], loadName[64];
    STImage *sourceImage, *targetImage;
    parseConfigFile(configFile.c_str(),
                    sourceName, targetName,
                    saveName, loadName,
                    &sourceImage, &targetImage);
    delete sourceImage;
    delete targetImage;

    //
    // load the features from the saved features file
    //
    loadLineEditorFile(loadName, AddFeatureCallback,
                       sourceName, targetName,
                       &sourceImage, &targetImage);

    //
    // run the full morphing algorithm before going into the main loop to
    // display an image
    //

    // these weighting parameters (Beier & Nelly 1992) can be changed if desired
    const float a = 0.5f, b = 1.0f, p = 0.2f;

    GenerateMorphFrames(sourceImage, gSourceFeatures,
                        targetImage, gTargetFeatures,
                        a, b, p);


    //
    // display a test or debug image here if desired
    // (note: comment this out if you call DisplayImage from elsewhere)
    //
    //STImage *result = sourceImage;

    // use this to test your image blending
    //STImage *result = BlendImages(sourceImage, targetImage, 0.2f);

    // use this to test your field morph
    STImage *result = FieldMorph(sourceImage, gSourceFeatures, gTargetFeatures,
                                 0.5f, a, b, p);

    // use this to test your image morphing
    /*
    STImage *result = MorphImages(sourceImage, gSourceFeatures,
                                  targetImage, gTargetFeatures,
                                  0.5f, a, b, p);
    */
    
    DisplayImage(result);

    // enter the GLUT main loop
    glutMainLoop();

    return 0;
}

// --------------------------------------------------------------------------
