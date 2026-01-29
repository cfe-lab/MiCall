
FROM debian:12

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq

RUN apt-get install -y ruby ruby-dev
RUN apt-get install -y build-essential

RUN gem install bundler

WORKDIR /docs

COPY ./Gemfile* .

RUN bundle install

EXPOSE 80

COPY . .

ENTRYPOINT ["bundle", "exec", "jekyll", "serve", "--incremental", "--livereload", "--host", "0.0.0.0", "--port", "80"]
