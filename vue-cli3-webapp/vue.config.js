module.exports = {
    publicPath: process.env.NODE_ENV === 'production'
        ? '/themes/custom/physics/mie/' //deploy path in Drupal setup at physics.ifmo.ru
        : '/'
}
