#include "glwidget.h"

#define STB_IMAGE_IMPLEMENTATION
#include "graphics/stb_image.h"
#include <QApplication>
#include <QKeyEvent>
#include <iostream>

#define SPEED 1.5
#define ROTATE_SPEED 0.0025

using namespace std;
using namespace Eigen;

GLWidget::GLWidget(QWidget *parent) :
    QOpenGLWidget(parent),
    m_sim(),
    m_camera(),
    m_defaultShader(),
    m_movementScaling(),
    // Movement
    m_deltaTimeProvider(),
    m_intervalTimer(),
    // Timing
    m_forward(),
    m_sideways(),
    m_vertical(),
    // Mouse handler stuff
    m_lastX(),
    m_lastY(),
    m_leftCapture(false),
    m_rightCapture(false),
    m_rightClickSelectMode(SelectMode::None),
    m_lastSelectedVertex(-1)
{
    std::cout << "glwidget constructor" << std::endl;

    // GLWidget needs all mouse move events, not just mouse drag events
    setMouseTracking(true);

    // Hide the cursor since this is a fullscreen app
    QApplication::setOverrideCursor(Qt::ArrowCursor);

    // GLWidget needs keyboard focus
    setFocusPolicy(Qt::StrongFocus);

    // Function tick() will be called once per interva
    connect(&m_intervalTimer, SIGNAL(timeout()), this, SLOT(tick()));
}

GLWidget::~GLWidget()
{
    if (m_defaultShader != nullptr) delete m_defaultShader;
}

// ================== Basic OpenGL Overrides

void GLWidget::initializeGL()
{
    // Initialize GL extension wrangler
    glewExperimental = GL_TRUE;
    GLenum err = glewInit();
    if (err != GLEW_OK) fprintf(stderr, "Error while initializing GLEW: %s\n", glewGetErrorString(err));
    fprintf(stdout, "Successfully initialized GLEW %s\n", glewGetString(GLEW_VERSION));

    // Set clear color to white
    glClearColor(1, 1, 1, 1);

    // Enable depth-testing and backface culling
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    // Initialize shaders
    m_defaultShader = new Shader(":resources/shaders/shader.vert", ":resources/shaders/shader.frag");
    m_skybox_shader = new Shader(":resources/shaders/skybox.vert", ":resources/shaders/skybox.frag"); // loads and compiles shaders

    setupSkybox();

    // Initialize ARAP, and get parameters needed to decide the camera position, etc
    Vector3f coeffMin, coeffMax;
    m_sim.init(coeffMin, coeffMax);

    Vector3f center = (coeffMax + coeffMin) / 2.0f;
    float extentLength  = (coeffMax - coeffMin).norm();

    // Scale all movement by this amount
    m_movementScaling = extentLength * 0.5;


    // Note for maintainers: Z-up
    float fovY = 120;
    float nearPlane = 0.0001f;
    float farPlane  = 100 * extentLength;

    // Initialize camera with a reasonable transform
    Eigen::Vector3f eye    = center - Eigen::Vector3f::UnitZ() * (extentLength-300) + Eigen::Vector3f::UnitY() * 300 + Eigen::Vector3f::UnitX() * 400;
    Eigen::Vector3f target = center;
    m_camera.lookAt(eye, target);
    m_camera.setOrbitPoint(target);
    m_camera.setPerspective(120, width() / static_cast<float>(height()), nearPlane, farPlane);

    m_deltaTimeProvider.start();
    m_intervalTimer.start(1000 / 60);
}

void GLWidget::setupSkybox() {

    print("hello");

    // from https://learnopengl.com/Advanced-OpenGL/Cubemaps
    std::vector<std::string> textures_faces = {
        "skyboxes/skybox5/right.jpg",
        "skyboxes/skybox5/left.jpg",
        "skyboxes/skybox5/top.jpg",
        "skyboxes/skybox5/bottom.jpg",
        "skyboxes/skybox5/front.jpg",
        "skyboxes/skybox5/back.jpg"
    };

    glGenTextures(1, &skyboxTexture);
    glBindTexture(GL_TEXTURE_CUBE_MAP, skyboxTexture);

    // from https://learnopengl.com/Advanced-OpenGL/Cubemaps
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    // from https://learnopengl.com/Advanced-OpenGL/Cubemaps
    for (unsigned int i = 0; i < 6; i++) {
        int width, height, nrChannels;
        unsigned char *data = stbi_load(textures_faces[i].c_str(), &width, &height, &nrChannels, 0);
        if (data) {
            // std::cout<<GL_TEXTURE_CUBE_MAP_POSITIVE_X + i<<std::endl;
            stbi_set_flip_vertically_on_load(false);
            glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X + i, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data);
            stbi_image_free(data);
        } else {
            std::cout << "Cubemap texture failed to load at path: " << textures_faces[i] << std::endl;
            stbi_image_free(data);
        }
    }

    //    std::cout<<GL_TEXTURE_CUBE_MAP_POSITIVE_X<<std::endl;
    //    std::cout<<GL_TEXTURE_CUBE_MAP_NEGATIVE_X<<std::endl;
    //    std::cout<<GL_TEXTURE_CUBE_MAP_POSITIVE_Y<<std::endl;
    //    std::cout<<GL_TEXTURE_CUBE_MAP_NEGATIVE_Y<<std::endl;
    //    std::cout<<GL_TEXTURE_CUBE_MAP_POSITIVE_Z<<std::endl;
    //    std::cout<<GL_TEXTURE_CUBE_MAP_NEGATIVE_Z<<std::endl;

    glBindTexture(GL_TEXTURE_CUBE_MAP, 0);

    GLfloat skyboxVertices[] = {

        1.f, 1.f, 1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f, -1.f, -1.f,  1.f, -1.f, 1.f, 1.f, 1.f, 1.f, 1.f,
        -1.f, 1.f, -1.f, -1.f, -1.f, -1.f, 1.f, -1.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, -1.f,
        1.f, 1.f, 1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f, -1.f, 1.f, 1.f, -1.f, 1.f, 1.f, 1.f,
        -1.f, -1.f, 1.f, 1.f, -1.f, 1.f, -1.f, -1.f, -1.f, -1.f, -1.f, -1.f, 1.f, -1.f, 1.f, 1.f, -1.f, -1.f,
        -1.f, 1.f, 1.f, -1.f, -1.f, 1.f, -1.f, -1.f, -1.f, -1.f, -1.f, -1.f, -1.f, 1.f, -1.f, -1.f, 1.f, 1.f,
        1.f, 1.f, 1.f, 1.f, 1.f, -1.f, 1.f, -1.f, -1.f, 1.f, -1.f, -1.f, 1.f, -1.f, 1.f, 1.f, 1.f, 1.f

    };

    glGenVertexArrays(1, &skyboxVAO);
    glGenBuffers(1, &skyboxVBO);

    glBindVertexArray(skyboxVAO);
    glBindBuffer(GL_ARRAY_BUFFER, skyboxVBO);
    glBufferData(GL_ARRAY_BUFFER, 108 * sizeof(GLfloat), &skyboxVertices, GL_STATIC_DRAW);

    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), (void*)0);

    glBindBuffer(GL_ARRAY_BUFFER, 0); // unbinds the buffer
    glBindVertexArray(0); // unbinds the vao
}

void GLWidget::paintGL()
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // ---------- SKYBOX ----------
    glDepthMask(GL_FALSE);
    glUseProgram(m_skybox_shader->getGLuint());

    glUniform1i(glGetUniformLocation(m_skybox_shader->getGLuint(), "u_skybox"), 0);

    // written by GPT to take top left 3x3 of 4x4 view matrix
    Eigen::Matrix4f tempMatrix = Eigen::Matrix4f::Identity();
    tempMatrix.block<3,3>(0,0) = m_camera.getView().block<3,3>(0,0);
    glUniformMatrix4fv(glGetUniformLocation(m_skybox_shader->getGLuint(), "u_viewMatrix"), 1, GL_FALSE, tempMatrix.data());

    glUniformMatrix4fv(glGetUniformLocation(m_skybox_shader->getGLuint(), "u_projectionMatrix"), 1, GL_FALSE, m_camera.getProjection().data());

    glBindVertexArray(skyboxVAO);
    glBindTexture(GL_TEXTURE_CUBE_MAP, skyboxTexture);

    glDrawArrays(GL_TRIANGLES, 0, 36);
    glDepthMask(GL_TRUE);

    // ---------- REFLECT AND REFRACT ----------

    m_defaultShader->bind();
    glUseProgram(m_defaultShader->getGLuint()); // use the shader program for rendering
    glUniform1i(glGetUniformLocation(m_defaultShader->getGLuint(), "u_skybox"), 0);
    glUniform1f(glGetUniformLocation(m_defaultShader->getGLuint(), "u_reflection"), 1);
    glUniform1f(glGetUniformLocation(m_defaultShader->getGLuint(), "u_refraction"), 0);
    glUniform1f(glGetUniformLocation(m_defaultShader->getGLuint(), "u_materialRefractiveIndex"), 0);
    glUniform1i(glGetUniformLocation(m_defaultShader->getGLuint(), "drawmode"), 0);
    m_defaultShader->setUniform("proj", m_camera.getProjection());
    m_defaultShader->setUniform("view", m_camera.getView());
    m_sim.draw(m_defaultShader, GL_TRIANGLES);

    glClear(GL_DEPTH_BUFFER_BIT);

    glUniform1i(glGetUniformLocation(m_defaultShader->getGLuint(), "drawmode"), 1);
    m_sim.drawShapes(m_defaultShader, GL_TRIANGLES);
    m_defaultShader->unbind();

    glClear(GL_DEPTH_BUFFER_BIT);


}

void GLWidget::resizeGL(int w, int h)
{
    glViewport(0, 0, w, h);
    m_camera.setAspect(static_cast<float>(w) / h);
}

// ================== Event Listeners

Eigen::Vector3f GLWidget::transformToWorldRay(int x, int y)
{
    Eigen::Vector4f clipCoords = Eigen::Vector4f(
        (float(x) / width()) * 2.f - 1.f,
        1.f - (float(y) / height()) * 2.f,
        -1.f,
        1.f);

    Eigen::Vector4f transformed_coords = m_camera.getProjection().inverse() * clipCoords;
    transformed_coords = Eigen::Vector4f(transformed_coords.x(), transformed_coords.y(), -1.f, 0.f);
    transformed_coords = m_camera.getView().inverse() * transformed_coords;

    return Eigen::Vector3f(transformed_coords.x(), transformed_coords.y(), transformed_coords.z()).normalized();
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{
    // Get current mouse coordinates
    const int currX = event->position().x();
    const int currY = event->position().y();

    // Get closest vertex to ray
    const Vector3f ray = transformToWorldRay(currX, currY);

    // Switch on button
    switch (event->button()) {

    case Qt::MouseButton::LeftButton: {
        // Capture
        m_leftCapture = true;
        // Select this vertex
        break;
    }
    default: break;
    }

    // Set last mouse coordinates
    m_lastX = currX;
    m_lastY = currY;
}

void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    // Return if neither mouse button is currently held down
    if (!(m_leftCapture || m_rightCapture)) {
        return;
    }

    // Get current mouse coordinates
    const int currX = event->position().x();
    const int currY = event->position().y();

    // Find ray
    const Vector3f ray = transformToWorldRay(event->position().x(), event->position().y());
    Vector3f pos;



    // If the selected point is an anchor point

    // Rotate the camera
    const int deltaX = currX - m_lastX;
    const int deltaY = currY - m_lastY;
    if (deltaX != 0 || deltaY != 0) {
        m_camera.rotate(deltaY * ROTATE_SPEED, -deltaX * ROTATE_SPEED);
    }


    // Set last mouse coordinates
    m_lastX = currX;
    m_lastY = currY;
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{
    m_leftCapture = false;
    m_lastSelectedVertex = -1;

    m_rightCapture = false;
    m_rightClickSelectMode = SelectMode::None;
}

void GLWidget::wheelEvent(QWheelEvent *event)
{
    float zoom = 1 - event->pixelDelta().y() * 0.1f / 120.f;
    m_camera.zoom(zoom);
}

void GLWidget::keyPressEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat()) return;

    switch (event->key())
    {
    case Qt::Key_W: m_forward  += SPEED; break;
    case Qt::Key_S: m_forward  -= SPEED; break;
    case Qt::Key_A: m_sideways -= SPEED; break;
    case Qt::Key_D: m_sideways += SPEED; break;
    case Qt::Key_F: m_vertical -= SPEED; break;
    case Qt::Key_R: m_vertical += SPEED; break;
    case Qt::Key_C: m_camera.toggleIsOrbiting(); break;
    case Qt::Key_Up: config.driveForce = Vector3d(1, 0, 0); break;
    case Qt::Key_Down: config.driveForce = Vector3d(-1, 0, 0); break;
    case Qt::Key_Left: config.rotateLeft = true; break;
    case Qt::Key_Right: config.rotateRight = true; break;
    case Qt::Key_Escape: QApplication::quit();
    }
}

void GLWidget::keyReleaseEvent(QKeyEvent *event)
{
    if (event->isAutoRepeat()) return;

    switch (event->key())
    {
    case Qt::Key_W: m_forward  -= SPEED; break;
    case Qt::Key_S: m_forward  += SPEED; break;
    case Qt::Key_A: m_sideways += SPEED; break;
    case Qt::Key_D: m_sideways -= SPEED; break;
    case Qt::Key_F: m_vertical += SPEED; break;
    case Qt::Key_R: m_vertical -= SPEED; break;
    case Qt::Key_Up: config.driveForce = Vector3d(0, 0, 0); break;
    case Qt::Key_Down: config.driveForce = Vector3d(0, 0, 0); break;
    case Qt::Key_Left: config.rotateLeft = false; break;
    case Qt::Key_Right: config.rotateRight = false; break;
    }
}

// ================== Physics Tick

void GLWidget::tick()
{
    float deltaSeconds = m_deltaTimeProvider.restart() / 1000.f;

    // Move camera
    auto look = m_camera.getLook();
    look.y() = 0;
    look.normalize();
    Eigen::Vector3f perp(-look.z(), 0, look.x());
    Eigen::Vector3f moveVec = m_forward * look.normalized() + m_sideways * perp.normalized() + m_vertical * Eigen::Vector3f::UnitY();
    moveVec *= m_movementScaling;
    moveVec *= deltaSeconds;
    m_camera.move(moveVec);
    m_sim.update((double)deltaSeconds);

    // Flag this view for repainting (Qt will call paintGL() soon after)
    update();

}
